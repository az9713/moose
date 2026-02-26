---
name: docker-windows
description: Run Docker containers on Windows with MINGW/MSYS2/Git Bash. Prevents path mangling, volume mount failures, and flag conversion issues that silently break Docker commands.
triggers:
  - docker
  - container
  - MOOSE
  - combined-opt
  - Docker Desktop
---
# Docker on Windows (MINGW/MSYS2/Git Bash)

## The Problem

Running Docker from Git Bash / MINGW on Windows silently breaks commands in three ways:

1. **Path mangling**: MINGW converts anything that looks like a Unix path to a Windows path.
   `-w /work` becomes `-w C:/Program Files/Git/work` — Docker gets a nonsensical Windows path.

2. **Flag conversion**: Single-letter flags like `-i` (for MOOSE input files) get eaten by MINGW's
   path-conversion heuristic when followed by a filename. `-i input.i` may be silently mangled.

3. **Volume mount paths**: Forward slashes in volume source paths can trigger path conversion.
   `-v /c/Users/...` may be rewritten to `-v C:/Program Files/Git/c/Users/...`.

These failures are **silent** — Docker may start, mount an empty volume, or pass garbled flags, and
the command appears to "work" while producing no output or wrong results.

## The Solution

**Always prefix Docker commands with `MSYS_NO_PATHCONV=1`** to disable all MINGW path conversion.

### Canonical Docker Run Pattern

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/USER/path/to/host/dir:/container/mount" \
  -w /container/mount/subdir \
  --entrypoint /bin/bash \
  IMAGE_NAME:TAG \
  -c 'COMMAND_TO_RUN 2>&1 | tail -N'
```

### Key Rules

| Rule | Why |
|------|-----|
| `MSYS_NO_PATHCONV=1` | Disables all MINGW path mangling for the command |
| Use `C:/Users/...` (forward slashes) for the host volume path | Docker on Windows accepts forward slashes; avoids escaping issues |
| Use `--entrypoint /bin/bash` with `-c '...'` | Wraps the entire command in a shell, preventing flag/path mangling inside the container |
| Pipe through `tail -N` | Captures only the last N lines of output (MOOSE is verbose) |
| Use `--rm` | Auto-removes the container after exit |
| Quote the `-c` argument with **single quotes** | Prevents the host shell from expanding variables or paths inside the container command |

### What NOT to Do

```bash
# BAD: No MSYS_NO_PATHCONV — paths will be mangled
docker run --rm -v "C:/mydir:/work" -w /work image:tag command

# BAD: Using -i flag directly — MINGW may eat it
docker run --rm -v "..." -w /work image:tag combined-opt -i input.i

# BAD: Using double quotes around -c — host shell expands paths
docker run --rm --entrypoint /bin/bash image:tag -c "cd /work && run"

# BAD: Using NUL instead of /dev/null inside container
# (The container is Linux even though the host is Windows)
```

## MOOSE-Specific Pattern

The MOOSE Docker image `idaholab/moose:latest` has `combined-opt` at `/opt/moose/bin/combined-opt`.

### Running a MOOSE Input File

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case01-1d-steady-diffusion \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case01_diffusion_1d.i 2>&1 | tail -30'
```

### Running Multiple Cases in Sequence

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '
    for dir in case14-thermoelasticity case15-lid-driven-cavity; do
      echo "=== $dir ==="
      cd /work/$dir
      ifile=$(ls *.i | head -1)
      /opt/moose/bin/combined-opt -i $ifile 2>&1 | tail -5
      echo
    done
  '
```

### MOOSE-Specific Gotchas

| Issue | Fix |
|-------|-----|
| `sh: mpicxx: command not found` / JIT failures | Add `disable_fpoptimizer = true` and `enable_jit = false` to `DerivativeParsedMaterial` blocks |
| Output files not appearing on host | Ensure the volume mount covers the working directory; Docker writes to the mounted path |
| `interval` parameter unused error | Use `time_step_interval` instead (MOOSE renamed this parameter) |
| PorousFlow missing material properties | `PorousFlowBasicTHM` does NOT auto-create `PorousFlowConstantBiotModulus` or `PorousFlowConstantThermalExpansionCoefficient` — add them as explicit `[Materials]` blocks |

## Debugging Checklist

If a Docker command produces no output or wrong results on Windows:

1. **Check `MSYS_NO_PATHCONV=1`** — is it set before the `docker` command?
2. **Check volume mount** — run `ls` inside the container to verify files are visible:
   ```bash
   MSYS_NO_PATHCONV=1 docker run --rm \
     -v "C:/Users/simon/mydir:/work" \
     --entrypoint /bin/bash image:tag \
     -c 'ls -la /work'
   ```
3. **Check Docker Desktop is running** — `docker info` should succeed
4. **Check entrypoint pattern** — use `--entrypoint /bin/bash` + `-c '...'` to avoid flag mangling
5. **Check output persistence** — files written inside the container only persist if they're in the mounted volume path

## General Docker Commands (Windows-Safe)

```bash
# Check Docker is running
docker info > /dev/null 2>&1 && echo "Docker OK" || echo "Start Docker Desktop"

# Pull an image
MSYS_NO_PATHCONV=1 docker pull idaholab/moose:latest

# Interactive shell (for debugging)
MSYS_NO_PATHCONV=1 docker run --rm -it \
  -v "C:/Users/simon/mydir:/work" \
  -w /work \
  idaholab/moose:latest \
  /bin/bash

# Run with resource limits (CPU-only machine)
MSYS_NO_PATHCONV=1 docker run --rm --cpus=4 --memory=8g \
  -v "C:/Users/simon/mydir:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c 'command here'
```
