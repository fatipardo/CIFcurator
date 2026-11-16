# Running the CIFcurator server on a cloud server (test / demo)

Everything needed to put the CIFcurator server (`../server/`) on a small rented Linux server, with
HTTPS, so that other people can use it from their web browser:

| File | What it is |
| --- | --- |
| `Dockerfile` | The image: Python, Biopython, numpy, pandas, MMseqs2 and MUSCLE (from `environment-server.yml`), and the CIFcurator code |
| `environment-server.yml` | The conda packages of the image (same versions as `../environment.yml`) |
| `docker-compose.yml` | Runs the server and Caddy (HTTPS) together; settings come from `.env` |
| `Caddyfile` | Caddy gets the HTTPS certificate by itself and passes the requests to the server |
| `.env.example` | Settings: address, limits per job, extra server options. `setup.sh` copies it to `.env` |
| `setup.sh` | Run once on a new server: installs Docker, firewall, security updates, swap, then starts everything |
| `update.sh` | Run after changing the code: rebuilds and restarts (jobs are kept) |
| `Dockerfile.dockerignore` | What is not sent to the image build (jobs, PDB sequence file, test cache, notebooks) |

The code itself is not in these files: changing the code only needs `update.sh`. Edit these files
only for a new dependency (a line in `environment-server.yml`) or a new server option (`.env`).

## 1. Rent a server

Any provider works. For example Hetzner Cloud (cheapest) or DigitalOcean:

- **Image:** Ubuntu 24.04
- **Size:** 4-8 vCPU, **16 GB RAM**, 160 GB disk (e.g. Hetzner CX42/CPX41, about 20-30 USD a month;
  DigitalOcean 4 vCPU / 16 GB, about 85 USD a month). A job can need 9 GB or more, and about
  1 GB of disk (kept 7 days). Billed by the hour: delete the server after the
  demo and you only pay for the days it existed.
- **SSH key:** add yours when creating the server (no password login).

The provider shows the server's IP address, e.g. `203.0.113.7`. Then:

```bash
ssh root@203.0.113.7
```

## 2. Put the code on the server

Either clone it from GitHub (if the server code is on a branch there):

```bash
git clone -b <branch> https://github.com/fatipardo/CIFcurator.git /opt/CIFcurator
```

or copy your checkout from your computer (from the folder that holds `server/` and `deploy/`):

```bash
rsync -a --exclude server/server_data --exclude tests/.cache --exclude '__pycache__' ./ root@203.0.113.7:/opt/CIFcurator/
```

## 3. Start it

```bash
cd /opt/CIFcurator/deploy
sudo bash setup.sh
```

It takes about 5-10 minutes the first time (installing the conda packages). At the end it prints the
address, e.g. **https://203-0-113-7.sslip.io**. sslip.io turns the IP address into a name, so the
server gets a real HTTPS certificate without buying a domain. The first job also downloads the PDB
sequence file (about 400 MB, once a week).

## Settings (`.env`)

| Setting | Default | Meaning |
| --- | --- | --- |
| `SITE_ADDRESS`, `PUBLIC_URL` | from the IP address | the server's address (a domain name later, e.g. `cifcurator.stanford.edu`) |
| `RETENTION_DAYS` | 7 | jobs are deleted this many days after they end |
| `MAX_ENTRIES` | 400 | a keyword that matches more PDB entries is refused before anything is downloaded |
| `MAX_CSM` | 200 | at most this many computed models (AlphaFold DB) per job |
| `JOB_TIMEOUT_HOURS` | 2 | a job still running after this many hours is stopped |
| `MEM_LIMIT` | 13g | memory of the CIFcurator container: a runaway job is killed, not the server |
| `EXTRA_ARGS` | (none) | other server options, e.g. `--show-recent-jobs` (list of everybody's jobs: fine among colleagues) or the email options `--smtp-host ... --mail-from ...` |

After editing `.env`: `docker compose up -d`.

## After changing the code

Get the new code on the server (`git pull` in `/opt/CIFcurator`, or the same `rsync` as above), then:

```bash
cd /opt/CIFcurator/deploy
bash update.sh
```

A job that is queued or running when the server restarts fails, so `update.sh` stops if there is one
(`bash update.sh --force` to restart anyway). Jobs, results and the PDB sequence file are kept (Docker
volume `cifcurator-data`). A new Python package or tool: add it to `environment-server.yml`, then
`update.sh`.

## Everyday commands (in `deploy/`)

```bash
docker compose ps                    # is it running?
docker compose logs -f cifcurator    # server log (jobs, emails)
docker compose logs -f caddy         # HTTPS / certificate problems
docker compose restart cifcurator
docker compose down                  # stop (jobs and certificate kept)
df -h /; free -h                     # disk and memory
```

The jobs are in the volume: `docker compose exec cifcurator ls /data/jobs`.

## When the demo is over

Delete the server in the provider's web page (this deletes the jobs too). Nothing else to clean up.

## Notes

- Docker opens ports 80 and 443 itself; the firewall (`ufw`) closes everything else. Port 8000 of the
  CIFcurator server is not reachable from outside: only Caddy talks to it.
- The server runs as an unprivileged user inside the container. Security updates of Ubuntu are
  installed automatically (`unattended-upgrades`). For a newer Caddy, change its version in
  `docker-compose.yml`, then `docker compose pull caddy && docker compose up -d`.
- Before a public (non-demo) server: per-person limits on queued jobs, a privacy/terms note, and,
  on a Stanford address, Stanford's security review.
