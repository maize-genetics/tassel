"""Build-time variables for the MkDocs site.

`build.gradle.kts` is the single source of truth for the released version, so
pages that quote the current release (Maven coordinates, standalone archive
names) use `{{ version }}` rather than a literal that has to be hand-edited
every release.

The Gradle file is parsed textually so the docs build needs no JDK: both
docs_ci.yml and deploy_project_site.yml run on Python-only runners.
"""

import re
from pathlib import Path

GRADLE_BUILD_FILE = Path(__file__).resolve().parent.parent / "build.gradle.kts"

_VERSION = re.compile(r'^version\s*=\s*"([^"]+)"', re.MULTILINE)
_VERSION_DATE = re.compile(r'^val\s+versionDate\s*=\s*"([^"]+)"', re.MULTILINE)


def _extract(pattern, text, label):
    match = pattern.search(text)
    if match is None:
        raise ValueError(f"Could not read {label} from {GRADLE_BUILD_FILE}")
    return match.group(1)


def define_env(env):
    text = GRADLE_BUILD_FILE.read_text(encoding="utf-8")
    env.variables["version"] = _extract(_VERSION, text, "the project version")
    env.variables["version_date"] = _extract(_VERSION_DATE, text, "versionDate")
