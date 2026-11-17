#!/bin/bash
# Run once on a fresh Ubuntu 24.04 server, from this folder:  sudo bash setup.sh
# Installs Docker (from Ubuntu's own packages; buildx so that Dockerfile.dockerignore is used),
# opens only SSH, HTTP and HTTPS, turns on automatic security updates, adds swap, writes .env (if
# missing) and starts CIFcurator + Caddy.
set -euo pipefail
cd "$(dirname "$0")"
[ "$(id -u)" = 0 ] || { echo "Run as root: sudo bash setup.sh"; exit 1; }

apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y docker.io docker-compose-v2 docker-buildx ufw unattended-upgrades git
systemctl enable --now docker
dpkg-reconfigure -f noninteractive unattended-upgrades

# Firewall: SSH, HTTP (certificate check + redirect) and HTTPS only
ufw allow OpenSSH
ufw allow 80/tcp
ufw allow 443/tcp
ufw --force enable

# 4 GB of swap, so that a large job slows down instead of being killed when memory runs out
if ! swapon --show | grep -q /swapfile; then
    fallocate -l 4G /swapfile && chmod 600 /swapfile && mkswap /swapfile && swapon /swapfile
    grep -q '^/swapfile' /etc/fstab || echo '/swapfile none swap sw 0 0' >> /etc/fstab
fi

# Settings: the address is made from the server's public IP address (203.0.113.7 -> 203-0-113-7.sslip.io)
if [ ! -f .env ]; then
    ip=$(ip -4 route get 1.1.1.1 | awk '{for (i = 1; i < NF; i++) if ($i == "src") print $(i + 1)}')
    name="${ip//./-}.sslip.io"
    sed -e "s|^SITE_ADDRESS=.*|SITE_ADDRESS=$name|" -e "s|^PUBLIC_URL=.*|PUBLIC_URL=https://$name|" .env.example > .env
    echo "Wrote .env for https://$name (edit it to change the limits or the address)"
fi
if [ -d ../.git ]; then
    sed -i "s|^CIFCURATOR_COMMIT=.*|CIFCURATOR_COMMIT=$(git -C .. describe --always --dirty)|" .env
fi

docker compose up -d --build
echo
echo "CIFcurator is starting at $(grep '^PUBLIC_URL=' .env | cut -d= -f2-)"
echo "The HTTPS certificate takes a minute the first time. Logs: docker compose logs -f"
