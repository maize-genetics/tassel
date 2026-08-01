#!/usr/bin/env bash
#
# Assembles and archives the TASSEL standalone distribution.
#
# Usage: build-standalone.sh <archive-name>
#
# <archive-name> is both the archive filename stem and the single top-level
# directory inside each archive, so `tassel-5-standalone-v5.2.98` produces
# dist/tassel-5-standalone-v5.2.98.{zip,tar.gz}, each unpacking into one
# tassel-5-standalone-v5.2.98/ directory. Releases through 5.2.97 archived the
# staging directory's *contents*, which dumped the jar, the launchers, and lib/
# straight into whichever directory the user happened to be standing in.
#
# Both the release and the nightly workflow call this script so the bundle is
# assembled exactly once anywhere. When the two workflows each kept their own
# copy of these steps the lists drifted, and 5.2.97 shipped a
# run_change_history.pl whose classpath still pointed into dist/.
#
# Expects `./gradlew build` to have already produced build/libs/sTASSEL.jar and
# build/libs/lib/.

set -euo pipefail

NAME="${1:?usage: build-standalone.sh <archive-name>}"
STAGE="dist/${NAME}"

for required in build/libs/sTASSEL.jar build/libs/lib; do
    if [ ! -e "$required" ]; then
        echo "FAIL: ${required} not found; run ./gradlew build first" >&2
        exit 1
    fi
done

# --- Stage the bundle -------------------------------------------------------------

rm -rf "$STAGE"
mkdir -p "$STAGE"
cp -r build/libs/lib "${STAGE}/"
cp build/libs/sTASSEL.jar "${STAGE}/"
cp scripts/*.pl scripts/*.bat "${STAGE}/"

# --- Point the launchers at the flattened layout ----------------------------------

# The in-repo scripts expect the Gradle tree, where the jar sits under dist/. In
# the bundle it sits beside the launchers. Globbing rather than listing filenames
# keeps a newly added launcher from being silently left unpatched.
sed -i.bak 's|\$top/dist/sTASSEL\.jar|$top/sTASSEL.jar|' "$STAGE"/*.pl
sed -i.bak 's|%TOP%\\dist\\sTASSEL\.jar|%TOP%\\sTASSEL.jar|' "$STAGE"/*.bat
rm -f "$STAGE"/*.bak

if grep -l 'dist[/\]sTASSEL\.jar' "$STAGE"/*.pl "$STAGE"/*.bat; then
    echo "FAIL: the launchers above still point at dist/sTASSEL.jar" >&2
    exit 1
fi

# --- Archive the directory, not its contents --------------------------------------

cd dist
rm -f "${NAME}.zip" "${NAME}.tar.gz"
zip -qr "${NAME}.zip" "$NAME"
tar -czf "${NAME}.tar.gz" "$NAME"

echo "Built ${NAME}.zip and ${NAME}.tar.gz in dist/"
