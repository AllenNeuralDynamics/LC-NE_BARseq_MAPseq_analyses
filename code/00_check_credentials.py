"""Start-of-run preflight: report Code Ocean API credential status up front.

08_generate_metadata.py needs the Code Ocean REST API (via the "Code Ocean API
Credentials" Secret) to write processing.json. This runs first and actually
calls the API so the credential status is known in the first ~30 seconds rather
than after the ~2-hour pipeline.

It never aborts the run: the analysis itself doesn't need the API, so a
credential-less reproducer can still run everything (processing.json is simply
skipped). If you intend to publish this run as a provenance-tracked asset and
see the warning below, Ctrl-C and attach the Secret before re-running.
"""

import sys

from co_provenance import fetch_co_provenance


def main() -> None:
    try:
        capsule_url, version = fetch_co_provenance()
    except RuntimeError as e:
        print(
            "\n"
            "  WARNING: Code Ocean API credentials not available --\n"
            f"    {e}\n"
            "    processing.json will be SKIPPED at the end of the run.\n"
            "    data_description.json will still be written.\n"
            "    Ctrl-C now and attach the 'Code Ocean API Credentials' Secret\n"
            "    if you intend to publish this run as a provenance-tracked asset.\n",
            file=sys.stderr,
        )
        return
    print(f"Preflight OK: credentials valid (capsule {capsule_url}, version {version}).")


if __name__ == "__main__":
    main()
