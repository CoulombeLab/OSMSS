#!/usr/bin/env bash
# Run this file to install OSMSS on Linux:  bash "Install OSMSS (Linux).sh"
cd "$(dirname "$0")" || exit 1
exec bash installer/unix/install-osmss.sh
