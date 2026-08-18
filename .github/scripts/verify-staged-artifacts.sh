#!/usr/bin/env bash
#
# Verifies the artifacts staged in build/staging-deploy before they are deployed to
# Maven Central. Maven Central is immutable, so a bad upload can only be corrected by
# publishing a new version.
#
# Release 5.2.97 shipped three byte-identical JARs (the main JAR published under the
# main, sources, and javadoc coordinates) and no fat JAR at all. Every check below
# exists to make that class of mistake fail the release instead of reaching users.
#
# Usage: verify-staged-artifacts.sh <version> [staging-dir]

set -euo pipefail

VERSION="${1:?usage: verify-staged-artifacts.sh <version> [staging-dir]}"
STAGING_ROOT="${2:-build/staging-deploy}"
DIR="${STAGING_ROOT}/net/maizegenetics/tassel/${VERSION}"

failures=0

fail() {
    echo "FAIL: $*" >&2
    failures=$((failures + 1))
}

pass() {
    echo "ok: $*"
}

if [ ! -d "$DIR" ]; then
    echo "FAIL: staging directory not found: $DIR" >&2
    exit 1
fi

echo "Verifying staged artifacts in $DIR"
echo

# --- Required files, each with a detached GPG signature ---------------------------

MAIN_JAR="tassel-${VERSION}.jar"
SOURCES_JAR="tassel-${VERSION}-sources.jar"
JAVADOC_JAR="tassel-${VERSION}-javadoc.jar"
FAT_JAR="tassel-${VERSION}-jar-with-dependencies.jar"
POM="tassel-${VERSION}.pom"

for artifact in "$MAIN_JAR" "$SOURCES_JAR" "$JAVADOC_JAR" "$FAT_JAR" "$POM"; do
    if [ ! -s "${DIR}/${artifact}" ]; then
        fail "missing or empty artifact: ${artifact}"
        continue
    fi
    if [ ! -s "${DIR}/${artifact}.asc" ]; then
        fail "missing GPG signature: ${artifact}.asc"
        continue
    fi
    pass "${artifact} present and signed"
done

# Nothing below can run meaningfully without the archives themselves.
if [ "$failures" -ne 0 ]; then
    echo
    echo "${failures} check(s) failed."
    exit 1
fi

# --- The four JARs must be distinct -----------------------------------------------

if command -v sha256sum >/dev/null 2>&1; then
    sha256() { sha256sum "$@"; }
else
    sha256() { shasum -a 256 "$@"; }
fi

distinct=$(cd "$DIR" && sha256 "$MAIN_JAR" "$SOURCES_JAR" "$JAVADOC_JAR" "$FAT_JAR" |
    awk '{print $1}' | sort -u | grep -c '')
if [ "$distinct" -ne 4 ]; then
    fail "expected 4 distinct JAR checksums, found ${distinct} (are the Jar tasks sharing an output file?)"
else
    pass "main, sources, javadoc, and fat JARs are all distinct"
fi

# --- Each JAR must contain what its classifier claims ------------------------------

# Listings are written to files first: `grep -q` closes the pipe on its first match,
# which would trip `pipefail` on the upstream `unzip`.
listing_dir=$(mktemp -d)
trap 'rm -rf "$listing_dir"' EXIT

listing() {
    local jar="$1"
    local out="${listing_dir}/${jar}.txt"
    [ -f "$out" ] || unzip -l "${DIR}/${jar}" > "$out"
    echo "$out"
}

if grep -q '\.class$' "$(listing "$SOURCES_JAR")"; then
    fail "${SOURCES_JAR} contains compiled classes"
elif ! grep -qE '\.(java|kt)$' "$(listing "$SOURCES_JAR")"; then
    fail "${SOURCES_JAR} contains no Java or Kotlin sources"
else
    pass "${SOURCES_JAR} contains sources only"
fi

if ! grep -q 'index\.html' "$(listing "$JAVADOC_JAR")"; then
    fail "${JAVADOC_JAR} has no index.html (did Dokka generation run?)"
else
    pass "${JAVADOC_JAR} contains generated documentation"
fi

if ! grep -q 'net/maizegenetics/tassel/TASSELMainApp\.class' "$(listing "$MAIN_JAR")"; then
    fail "${MAIN_JAR} does not contain the TASSEL classes"
else
    pass "${MAIN_JAR} contains the TASSEL classes"
fi

# The fat JAR must carry its dependencies, not just TASSEL's own classes.
if ! grep -q 'org/apache/logging/log4j/' "$(listing "$FAT_JAR")"; then
    fail "${FAT_JAR} does not bundle its dependencies"
else
    pass "${FAT_JAR} bundles its dependencies"
fi

echo
if [ "$failures" -ne 0 ]; then
    echo "${failures} check(s) failed."
    exit 1
fi

echo "All staged artifacts look correct."
