#!/bin/bash
# Double-click this file on a Mac to install OSMSS (it opens in Terminal).
# If macOS blocks it, right-click it, choose Open, then Open again.
cd "$(dirname "$0")" || exit 1
exec /bin/bash installer/unix/install-osmss.sh
