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
      The asset name is built from ASSET_BASE_NAME + creation time, per the
      AIND create_processing_metadata guidance for derived assets.
  - processing.json : describes this analysis run (capsule URL + release
      version pulled from the Code Ocean REST API at runtime).

All schema construction goes through aind-data-schema's Pydantic models, so any
schema violation raises and aborts the run.
"""

from datetime import datetime, timezone
from pathlib import Path

from aind_data_schema.components.identifiers import Code, DataAsset
from aind_data_schema.core.data_description import DataDescription
from aind_data_schema.core.processing import DataProcess, ProcessStage, Processing
from aind_data_schema_models.data_name_patterns import DataLevel, build_data_name
from aind_data_schema_models.process_names import ProcessName

from co_provenance import fetch_co_provenance

RESULTS_DIR = Path("/results")
DATA_DIR = Path("/data")
EXPERIMENTERS = ["Polina Kosillo"]
ASSET_BASE_NAME = "BARseq-MAPseq-LC-NE-combined"
AIND_OPEN_DATA_BUCKET = "s3://aind-open-data"


def find_input_descriptions() -> list[tuple[Path, DataDescription]]:
    """Discover attached input assets and load their data_description.json.

    Returns (asset_dir, DataDescription) for every subdirectory of /data/ that
    contains a data_description.json. The attached assets are the four analysis
    inputs (2 BARseq + 2 MAPseq), so this finds exactly those.
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


def write_data_description(inputs: list[tuple[Path, DataDescription]], creation_time: datetime) -> None:
    """Build and write the DERIVED data_description.json for the analysis output.

    Constructed directly (rather than via DataDescription.from_data_description)
    so the asset name comes from ASSET_BASE_NAME instead of a single input's
    lineage -- this asset derives from four inputs across two subjects and two
    modalities, so an input-derived name would be misleading. Institution,
    funding, investigators and project_name are inherited from the first input.
    """
    descriptions = [dd for _, dd in inputs]
    source_names = sorted(dd.name for dd in descriptions)
    base_dd = descriptions[0]

    derived_dd = DataDescription(
        name=build_data_name(ASSET_BASE_NAME, creation_datetime=creation_time),
        creation_time=creation_time,
        institution=base_dd.institution,
        funding_source=base_dd.funding_source,
        investigators=base_dd.investigators,
        project_name=base_dd.project_name,
        modalities=dedupe_modalities(descriptions),
        data_level=DataLevel.DERIVED,
        subject_id=None,  # spans two subjects; both captured in source_data
        source_data=source_names,
        data_summary=(
            "Combined BARseq/MAPseq analysis of LC-NE neurons across subjects "
            "780345 and 780346 (manuscript Figure S5)."
        ),
    )
    derived_dd.write_standard_file(output_directory=RESULTS_DIR)
    print(f"  wrote data_description.json (name: {derived_dd.name})")
    print(f"    source_data: {source_names}")


def write_processing(inputs: list[tuple[Path, DataDescription]], start_time: datetime) -> None:
    """Build and write processing.json describing this analysis run.

    Skipped with a warning if the Code Ocean API credentials are not available
    (provenance fields can't be populated reliably without them); in that case
    only data_description.json is written. The start-of-run preflight
    (00_check_credentials.py) flags the same condition up front.
    """
    try:
        capsule_url, version = fetch_co_provenance()
    except RuntimeError as e:
        print(f"  WARNING: skipping processing.json -- {e}")
        return

    input_assets = [DataAsset(url=f"{AIND_OPEN_DATA_BUCKET}/{dd.name}") for _, dd in inputs]
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
        start_date_time=start_time,
        notes="BARseq normalization/clustering and BARseq-MAPseq projection analysis (Figure S5).",
    )
    processing = Processing(data_processes=[process])
    processing.write_standard_file(output_directory=RESULTS_DIR)
    print("  wrote processing.json")


def main() -> None:
    """Generate data_description.json and processing.json in /results/."""
    now = datetime.now(timezone.utc)
    inputs = find_input_descriptions()
    print(f"Found {len(inputs)} input asset(s) with data_description.json")
    for asset_dir, _ in inputs:
        print(f"  input: {asset_dir.name}")
    write_data_description(inputs, now)
    write_processing(inputs, now)
    print("Done.")


if __name__ == "__main__":
    main()
