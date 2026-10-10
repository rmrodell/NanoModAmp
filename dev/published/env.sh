#!/bin/bash
# Source this before running any dev/published preprocessing or counting script (no logic changes).
# Clean module environment: an interactive session's own modules make `ml py-cutadapt/1.18_py36` fail.
if ! type module >/dev/null 2>&1; then source /etc/profile >/dev/null 2>&1; fi
module --force purge >/dev/null 2>&1
module load devel math >/dev/null 2>&1

# Tool binaries (see versions.txt). Override if relocated.
export MINIMAP2="${MINIMAP2:-/oak/stanford/groups/nicolemm/rodell/minimap2/minimap2}"      # 2.28-r1221-dirty
export SEQTK="${SEQTK:-/home/users/rodell/seqtk/seqtk}"                                       # 1.5-r133
export UMICOLLAPSE="${UMICOLLAPSE:-/home/users/rodell/UMICollapse/umicollapse}"               # efeab35; java/11 (scripts load it)
# umi_tools 1.1.6 is a pip --user install under python/3.6.1 (~/.local/bin); the scripts load python/3.6.1.
case ":$PATH:" in *":$HOME/.local/bin:"*) ;; *) export PATH="$PATH:$HOME/.local/bin" ;; esac
