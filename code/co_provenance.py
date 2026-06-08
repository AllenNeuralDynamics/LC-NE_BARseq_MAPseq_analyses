"""Shared Code Ocean provenance lookup (no third-party dependencies).

Both the start-of-run credentials preflight (00_check_credentials.py) and the
metadata generator (08_generate_metadata.py) import this, so they validate and
use the exact same logic. urllib/json/base64 only -- deliberately free of
aind-data-schema so the preflight stays lightweight.
"""

import base64
import json
import os
import urllib.error
import urllib.request

CO_API_BASE = "https://codeocean.allenneuraldynamics.org/api/v1"
CO_WEB_BASE = "https://codeocean.allenneuraldynamics.org/capsule"


def fetch_co_provenance() -> tuple[str, str]:
    """Return (capsule_url, version) for the running Code Ocean capsule.

    Calls the Code Ocean REST API at runtime to look up the capsule's web URL
    (built from the slug) and the release version of this run. ``version`` is
    ``"from non-release editable capsule"`` when running an editable capsule;
    otherwise a string like ``"v1.0"``.

    Requires the "Code Ocean API Credentials" Secret to be attached to the
    capsule (Capsule Settings -> Credentials), which exposes API_KEY at runtime.
    Raises RuntimeError if any required env var is missing or the API call fails
    (e.g. the wrong/expired credentials are attached).
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
