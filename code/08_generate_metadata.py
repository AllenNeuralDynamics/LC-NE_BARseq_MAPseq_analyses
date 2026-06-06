"""Generate AIND-compliant derived-asset metadata for this analysis capsule.

Writes two files into /results/ so the run output can be published as a
DERIVED data asset in aind-open-data (a downstream capsule consumes it as
input, so it needs provenance):

  - data_description.json : DERIVED description of the combined analysis.
      Inherits institution / funding / investigators / project_name from one
      input asset's data_description.json, sets data_level=DERIVED, lists all
      input assets in source_data, and leaves subject_id unset (this asset
      spans two subjects -- 780345 and 780346 -- which the schema's single
      subject_id field can't represent; both are captured via source_data).
  - processing.json : describes this analysis run (capsule URL + release
      version pulled from the Code Ocean REST API at runtime).

Mirrors the metadata step in the sibling capsule LC-NE_BARseq_MAT-RDS_conversion.
All schema construction goes through aind-data-schema's Pydantic models, so any
schema violation raises and aborts the run.
"""

import base64
import json
import os
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from aind_data_schema.components.identifiers import Code, DataAsset
from aind_data_schema.core.data_description import DataDescription
from aind_data_schema.core.processing import DataProcess, ProcessStage, Processing
from aind_data_schema_models.process_names import ProcessName

RESULTS_DIR = Path("/results")
DATA_DIR = Path("/data")
EXPERIMENTERS = ["Polina Kosillo"]
PROCESS_NAME_LABEL = "analyzed"  # becomes the process token in the derived name

CO_API_BASE = "https://codeocean.allenneuraldynamics.org/api/v1"
CO_WEB_BASE = "https://codeocean.allenneuraldynamics.org/capsule"
AIND_OPEN_DATA_BUCKET = "s3://aind-open-data"


def fetch_co_provenance() -> tuple[str, str]:
    """Return (capsule_url, version) for the running Code Ocean capsule.

    Calls the Code Ocean REST API at runtime to look up the capsule's web URL
    (built from the slug) and the release version of this run. ``version`` is
    ``"from non-release editable capsule"`` when running an editable capsule;
    otherwise a string like ``"v1.0"``.

    Requires the "Code Ocean API Credentials" Secret to be attached to the
    capsule (Capsule Settings -> Credentials), which exposes API_KEY at
    runtime. Raises RuntimeError if any required env var is missing or the API
    call fails.
    """
    api_key = os.environ.get("API_KEY")
    capsule_id = os.environ.get("CO_CAPSULE_ID")
    computation_id = os.environ.get("CO_COMPUTATION_ID")
    if not api_key or not capsule_id or not computation_id:
        raise RuntimeError(
            "Missing Code Ocean env vars (API_KEY / CO_CAPSULE_ID / "
            "CO_COMPUTATION_ID). Attach the 'Code Ocean API Credentials' "
            "Secret to the capsule (Capsule Settings -> Credentials)."
        )

    auth = base64.b64encode(f"{api_key}:".encode()).decode()
    headers = {"Authorization": f"Basic {auth}"}

    def _get(path: str) -> dict:
        req = urllib.request.Request(f"{CO_API_BASE}{path}", headers=headers)
        with urllib.request.urlopen(req, timeout=10) as resp:
            return json.loads(resp.read())

    try:
        capsule = _get(f"/capsules/{capsule_id}")
        computation = _get(f"/computations/{computation_id}")
    except (urllib.error.URLError, json.JSONDecodeError) as e:
        raise RuntimeError(f"Code Ocean API call failed: {e}") from e

    capsule_url = f"{CO_WEB_BASE}/{capsule['slug']}/tree"
    if "version" in computation:
        version = f"v{computation['version']}.0"
    else:
        version = "from non-release editable capsule"
    return capsule_url, version


def find_input_descriptions() -> list[tuple[Path, DataDescription]]:
    """Discover attached input assets and load their data_description.json.

    Returns (asset_dir, DataDescription) for every subdirectory of /data/ that
    contains a data_description.json. After PR4 the only attached assets are
    the four analysis inputs (2 BARseq + 2 MAPseq), so this finds exactly those.
    """
    found = []
    for asset_dir in sorted(p for p in DATA_DIR.iterdir() if p.is_dir()):
        dd_path = asset_dir / "data_description.json"
        if not dd_path.exists():
            continue
        dd = DataDescription.model_validate_json(dd_path.read_text())
        found.append((asset_dir, dd))
    if not found:
        raise RuntimeError(f"No input assets with data_description.json under {DATA_DIR}")
    return found


def dedupe_modalities(descriptions: list[DataDescription]) -> list:
    """Union the modalities across all input assets, preserving order."""
    seen = set()
    union = []
    for dd in descriptions:
        for mod in dd.modalities:
            key = getattr(mod, "abbreviation", str(mod))
            if key not in seen:
                seen.add(key)
                union.append(mod)
    return union


def write_data_description(inputs: list[tuple[Path, DataDescription]]) -> None:
    """Build and write the DERIVED data_description.json for the analysis output."""
    descriptions = [dd for _, dd in inputs]
    source_names = sorted(dd.name for dd in descriptions)
    base_dd = descriptions[0]  # inherit institution/funding/investigators/project_name

    derived_dd = DataDescription.from_data_description(
        base_dd,
        process_name=PROCESS_NAME_LABEL,
        source_data=source_names,
        modalities=dedupe_modalities(descriptions),
        data_summary=(
            "Combined BARseq/MAPseq analysis of LC-NE neurons across subjects "
            "780345 and 780346 (manuscript Figure S5)."
        ),
    )
    # This asset spans two subjects; the single subject_id field can't hold both,
    # so clear it. Provenance for both subjects is captured in source_data.
    derived_dd = derived_dd.model_copy(update={"subject_id": None})
    derived_dd.write_standard_file(output_directory=RESULTS_DIR)
    print(f"  wrote data_description.json (source_data: {source_names})")


def write_processing(inputs: list[tuple[Path, DataDescription]]) -> None:
    """Build and write processing.json describing this analysis run.

    Skipped with a warning if the Code Ocean API credentials are not available
    (provenance fields can't be populated reliably without them).
    """
    try:
        capsule_url, version = fetch_co_provenance()
    except RuntimeError as e:
        print(f"  WARNING: skipping processing.json -- {e}")
        return

    input_assets = [
        DataAsset(url=f"{AIND_OPEN_DATA_BUCKET}/{dd.name}") for _, dd in inputs
    ]
    code = Code(
        url=capsule_url,
        name="LC-NE_BARseq_MAPseq_analyses",
        version=version,
        run_script=Path("code/run"),
        language="R",
        input_data=input_assets,
    )
    process = DataProcess(
        process_type=ProcessName.ANALYSIS,
        name="LC-NE BARseq/MAPseq analysis",
        stage=ProcessStage.ANALYSIS,
        code=code,
        experimenters=EXPERIMENTERS,
        start_date_time=datetime.now(timezone.utc),
        notes="BARseq normalization/clustering and BARseq-MAPseq projection analysis (Figure S5).",
    )
    processing = Processing(data_processes=[process])
    processing.write_standard_file(output_directory=RESULTS_DIR)
    print("  wrote processing.json")


def main() -> None:
    """Generate data_description.json and processing.json in /results/."""
    inputs = find_input_descriptions()
    print(f"Found {len(inputs)} input asset(s) with data_description.json")
    for asset_dir, _ in inputs:
        print(f"  input: {asset_dir.name}")
    write_data_description(inputs)
    write_processing(inputs)
    print("Done.")


if __name__ == "__main__":
    main()
