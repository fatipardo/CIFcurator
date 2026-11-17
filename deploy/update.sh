#!/bin/bash
# After the code changed (git pull, or files copied with rsync): rebuild the image and restart the
# server. Run from this folder:  bash update.sh   (sudo if your user is not in the docker group)
# Jobs, results and the PDB sequence file are kept. A job that is queued or running when the
# server restarts fails ("The server was stopped before this job finished"), so this script
# stops if there is one, unless you add --force.
set -euo pipefail
cd "$(dirname "$0")"
# python is only on the PATH through the image's entrypoint (micromamba), also for docker compose exec
busy=$(docker compose exec -T cifcurator /usr/local/bin/_entrypoint.sh python -c "import glob, json; print(sum(json.load(open(p))['state'] in ('queued', 'running') for p in glob.glob('/data/jobs/*/job.json')))" 2>/dev/null) || busy=0
if ! [[ "$busy" =~ ^[0-9]+$ ]]; then
    echo "Could not count the jobs of the server. If none is queued or running: bash update.sh --force"
    [ "${1:-}" = --force ] || exit 1
elif [ "$busy" != 0 ] && [ "${1:-}" != --force ]; then
    echo "$busy job(s) queued or running: they would fail. Try again later, or: bash update.sh --force"
    exit 1
fi
if [ -d ../.git ]; then
    sed -i "s|^CIFCURATOR_COMMIT=.*|CIFCURATOR_COMMIT=$(git -C .. describe --always --dirty)|" .env
fi
docker compose up -d --build
docker image prune -f > /dev/null
echo "Updated: $(grep '^CIFCURATOR_COMMIT=' .env | cut -d= -f2-)"
