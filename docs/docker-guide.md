# Running MOOSE with Docker: A Complete Guide

This guide is written for Windows users who have never used Docker before and want to run
MOOSE simulations without the pain of compiling a large scientific computing framework on
Windows. By the end of this guide you will be able to pull a pre-built MOOSE image, mount
your input files into the container, and produce Exodus output files on your Windows
filesystem — all with copy-paste commands.

No prior Linux, Docker, or MOOSE experience is assumed.

---

## Table of Contents

1. [Why MOOSE Needs Docker (Especially on Windows)](#1-why-moose-needs-docker-especially-on-windows)
2. [What is Docker?](#2-what-is-docker)
3. [Installing Docker Desktop on Windows](#3-installing-docker-desktop-on-windows)
4. [Pulling the MOOSE Image](#4-pulling-the-moose-image)
5. [Running MOOSE Simulations](#5-running-moose-simulations)
6. [Docker for Custom MOOSE Applications](#6-docker-for-custom-moose-applications)
7. [Troubleshooting](#7-troubleshooting)
8. [Docker Commands Reference](#8-docker-commands-reference)
9. [Comparison: Docker vs WSL2 vs Native Linux](#9-comparison-docker-vs-wsl2-vs-native-linux)

---

## 1. Why MOOSE Needs Docker (Especially on Windows)

### MOOSE is a Linux-Native Framework

MOOSE (Multiphysics Object-Oriented Simulation Environment) is a C++ scientific computing
framework developed at Idaho National Laboratory. It is designed to run on Linux and macOS
systems used in high-performance computing (HPC) environments. There is no native Windows
port, and building one from scratch is not practical for most users.

The reason comes down to MOOSE's dependencies. A working MOOSE installation requires all
of the following to be compiled and correctly linked together:

| Dependency | Purpose |
|------------|---------|
| GCC / Clang | C++ compiler with C++17 support |
| MPI (OpenMPI or MPICH) | Parallel communication between solver processes |
| PETSc | Nonlinear and linear solver library (the core numerical engine) |
| libMesh | Finite element library (mesh handling, shape functions, integration) |
| HDF5 | Binary output file format |
| METIS / ParMETIS | Parallel mesh partitioning |
| Python 3 | Test harness, documentation, utilities |
| NetCDF | Scientific data format (Exodus output) |

Building all of these from source on Linux takes 30-90 minutes and involves careful version
matching. On Windows, the situation is significantly worse: none of these libraries ship
official Windows binaries in the right configuration, and the MPI implementations available
on Windows (Microsoft MPI, Intel MPI) use different ABIs than the Linux versions MOOSE
expects. PETSc's Windows support is limited and not maintained in the MOOSE ecosystem.

### The Three Options for Windows Users

You have exactly three practical paths to running MOOSE on a Windows machine:

**Option 1 — A separate Linux machine or HPC cluster.**
The best performance. Used by researchers at national laboratories. Not accessible to
most students or engineers working on laptops.

**Option 2 — Windows Subsystem for Linux 2 (WSL2).**
Microsoft's Linux compatibility layer, now built into Windows 10/11. WSL2 runs a real
Linux kernel inside a lightweight virtual machine. MOOSE can be compiled and run in WSL2
with full MPI support. The downside is that you still need to compile MOOSE yourself
(30-90 minutes) and manage a Linux environment. Covered in the comparison table at the
end of this document.

**Option 3 — Docker.**
Docker packages a complete, pre-built MOOSE environment into a "container." You download
the pre-built image (no compilation) and run simulations immediately. Your input files
live on Windows and are made available to the container through a shared folder. Output
files written by MOOSE appear back on your Windows filesystem automatically.

This guide covers Option 3.

### The Architecture: What Runs Where

```
+--------------------------------------------------------------------+
|  Windows Host Machine                                              |
|                                                                    |
|  C:\Users\you\project\                                             |
|    input.i          (you write this)                               |
|    output_out.e     (MOOSE writes this back)                       |
|    output.csv       (MOOSE writes this back)                       |
|                                                                    |
|  +--------------------------------------------------------------+  |
|  |  Docker Container (Linux environment)                        |  |
|  |                                                              |  |
|  |  +------------+  +----------+  +-----------+  +----------+  |  |
|  |  |  MOOSE     |  |  PETSc   |  |  libMesh  |  |   MPI    |  |  |
|  |  |  framework |  | (solvers)|  | (elements)|  |(parallel)|  |  |
|  |  +------------+  +----------+  +-----------+  +----------+  |  |
|  |                                                              |  |
|  |  +------------+  +----------+  +-----------+  +----------+  |  |
|  |  |   HDF5     |  |  NetCDF  |  |  Python3  |  |  GCC/G++ |  |  |
|  |  | (file I/O) |  | (Exodus) |  |  (tools)  |  |(compiler)|  |  |
|  |  +------------+  +----------+  +-----------+  +----------+  |  |
|  |                                                              |  |
|  |  /work/  <===  shared volume  ===>  C:\Users\you\project\  |  |
|  |                                                              |  |
|  +--------------------------------------------------------------+  |
|                                                                    |
+--------------------------------------------------------------------+
```

The key insight from this diagram: your input files (`.i`) and output files (`.e`, `.csv`)
live permanently on your Windows filesystem. The container only borrows them temporarily
during the simulation run. When the container exits, the output files remain on Windows.
The container itself is stateless and is deleted after each run (using `--rm`).

---

## 2. What is Docker?

### Containers in Plain English

Imagine you want to share a recipe with a friend, but instead of a recipe card, you hand
them a fully equipped kitchen with every ingredient already measured out and every appliance
already pre-heated. They just press "go." That is what a Docker container does for
software.

A container is a lightweight, isolated Linux environment that includes:
- A minimal Linux operating system
- All installed software and libraries
- All configuration already done

Containers are not full virtual machines. They share the host computer's kernel (the core
of the operating system) but have their own isolated filesystem, network, and process
space. This makes them faster to start and less resource-intensive than a traditional VM
like VirtualBox or VMware.

### Images vs. Containers

These two words are used constantly in Docker and it is important to understand the
difference:

| Term | What it is | Analogy |
|------|-----------|---------|
| **Image** | A read-only snapshot of a complete environment, stored on disk | A DVD or ISO file |
| **Container** | A running instance created from an image | A movie playing from that DVD |

One image can spawn many containers simultaneously. Each container is independent — changes
in one do not affect others or the original image.

When you run `docker run idaholab/moose:latest`, Docker:
1. Looks for the `idaholab/moose` image locally
2. If not found, downloads it from Docker Hub
3. Creates a new container from that image
4. Runs your command inside the container
5. (With `--rm`) Deletes the container when done

### Docker Hub

Docker Hub (hub.docker.com) is a public registry where organizations publish pre-built
images. It works like an app store for containers. Idaho National Laboratory publishes
official MOOSE images there under the `idaholab` organization.

You do not need an account on Docker Hub to pull public images.

### The Docker Workflow

```
  Docker Hub                  Your Machine
  +---------------+           +----------------+   run    +---------------+
  |  idaholab/    |  docker   |  idaholab/     | ──────►  |  Container    |
  |  moose:latest |  pull     |  moose:latest  |          |  (running)    |
  |  (remote)     | ────────► |  (local copy)  |          |               |
  +---------------+           +----------------+          +-------+-------+
                                                                  |
                                                          /work/ directory
                                                          inside container
                                                                  |
                                                          ↕ shared volume
                                                                  |
                                                    C:\Users\you\simulations\
                                                    on your Windows filesystem
```

Once the image is pulled (downloaded once), subsequent `docker run` commands start
immediately without any network access needed.

---

## 3. Installing Docker Desktop on Windows

Docker Desktop is the official application that installs and manages Docker on Windows. It
provides a graphical dashboard, but you will mostly use it through the command line.

### 3.1 System Requirements

Before installing, verify your system meets these requirements:

**Minimum requirements:**
- Windows 10 version 2004 (build 19041) or later, or Windows 11
- 64-bit processor
- 4 GB RAM (8 GB recommended for MOOSE)
- BIOS-level virtualization enabled (see note below)

**Recommended for MOOSE:**
- 8 GB RAM or more
- 4+ CPU cores
- 20 GB free disk space (image is 2-4 GB; leave room for output files)

**Note on virtualization:** Docker requires hardware virtualization. On most modern
computers this is enabled by default. If Docker fails to install with a virtualization
error, you need to enter your computer's BIOS/UEFI settings and enable "Intel VT-x" or
"AMD-V." Search for your specific laptop/motherboard model plus "enable virtualization
BIOS" for instructions.

**Check your Windows version:**
Press `Win + R`, type `winver`, press Enter. The version number appears in the "Version"
line. You need 2004 or higher.

**Check if virtualization is enabled:**
Open Task Manager (`Ctrl + Shift + Esc`), click the "Performance" tab, select "CPU." Look
for "Virtualization: Enabled" in the right panel.

### 3.2 Enabling WSL2 (Recommended Backend)

Docker Desktop on Windows can use either Hyper-V or WSL2 as its backend. WSL2 is
strongly recommended because:
- Better file system performance when accessing Windows files
- Lower memory overhead
- Does not require Windows Pro/Enterprise (Hyper-V does)
- More actively maintained

To enable WSL2 before installing Docker:

1. Open PowerShell as Administrator (right-click Start, select "Windows PowerShell (Admin)")

2. Run these two commands, one at a time:
   ```powershell
   dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
   dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
   ```

3. Restart your computer.

4. After restart, open PowerShell again and run:
   ```powershell
   wsl --set-default-version 2
   ```

5. If you see a message about the WSL2 kernel update, download and install it from:
   `https://aka.ms/wsl2kernel`

6. Run step 4 again after installing the kernel update.

### 3.3 Downloading Docker Desktop

1. Open a web browser and navigate to: `https://www.docker.com/products/docker-desktop/`

2. Click "Download for Windows." The downloaded file will be named something like
   `Docker Desktop Installer.exe` (approximately 500 MB).

3. Do not run the installer yet — finish reading this section first.

### 3.4 Running the Installer

1. Double-click `Docker Desktop Installer.exe` to start the installer.

2. If prompted by Windows User Account Control, click "Yes."

3. On the "Configuration" screen:
   - Leave "Use WSL 2 instead of Hyper-V (recommended)" checked
   - Leave "Add shortcut to desktop" checked if desired
   - Click "Ok"

4. The installer will unpack files and configure the system. This takes 2-5 minutes.

5. When installation completes, click "Close and restart" if prompted, or "Close" and
   restart manually.

### 3.5 First Launch

After restarting:

1. Docker Desktop should start automatically. If not, find it in the Start menu or
   on your desktop and launch it.

2. The first launch involves additional setup. Accept the Docker Subscription Service
   Agreement (Docker Desktop is free for personal and educational use).

3. Docker Desktop will show a tutorial. You can skip it.

4. Wait for the Docker Desktop dashboard to show "Engine running" in the bottom left
   corner. This can take 30-60 seconds.

5. The Docker icon (a whale carrying containers) should appear in your system tray.

### 3.6 Verifying the Installation

Open a terminal. On Windows you can use:
- **Git Bash** (recommended for this guide — install from git-scm.com if you do not have it)
- **Command Prompt** (cmd.exe)
- **PowerShell**
- **Windows Terminal** (from the Microsoft Store)

Run the following commands:

**Check Docker version:**
```bash
docker --version
```

Expected output (version numbers will vary):
```
Docker version 25.0.3, build 4debf41
```

**Verify the Docker engine is running:**
```bash
docker run hello-world
```

Expected output:
```
Unable to find image 'hello-world:latest' locally
latest: Pulling from library/hello-world
...
Hello from Docker!
This message shows that your installation appears to be working correctly.
...
```

If you see this message, Docker is installed and working. The "Unable to find image
locally" message is normal — Docker downloaded the tiny `hello-world` image from Docker
Hub to verify the connection works.

### 3.7 Configuring Docker Desktop Resources

By default Docker Desktop may not allocate enough memory for large MOOSE simulations. To
adjust:

1. Open Docker Desktop
2. Click the gear icon (Settings) in the top right
3. Click "Resources" in the left sidebar
4. Set Memory to at least **4 GB** (8 GB recommended for complex problems)
5. Set CPUs to at least **2** (match the number of MPI processes you intend to use)
6. Click "Apply & Restart"

### 3.8 Enabling File Sharing (Important)

Docker Desktop must be allowed to share folders from your Windows filesystem with
containers. By default on WSL2, the entire Windows filesystem is accessible. If you use
Hyper-V backend, you may need to explicitly add folders:

1. Open Docker Desktop Settings
2. Click "Resources" then "File Sharing"
3. Add the drive or folder where your MOOSE input files will live (e.g., `C:\`)
4. Click "Apply & Restart"

### 3.9 Common Installation Issues

**Issue: "WSL2 installation is incomplete"**
Solution: Install the WSL2 Linux kernel update package from `https://aka.ms/wsl2kernel`
and restart.

**Issue: "Hardware assisted virtualization and data execution protection must be enabled"**
Solution: Reboot into BIOS/UEFI and enable VT-x (Intel) or AMD-V (AMD) in the CPU or
Advanced settings. The exact menu location varies by manufacturer.

**Issue: Docker Desktop shows "Engine stopped" and will not start**
Solution: Open PowerShell as Administrator and run:
```powershell
wsl --shutdown
```
Then restart Docker Desktop from the system tray.

**Issue: Installation fails with "HyperVisorLaunchType is off"**
Solution: Open Command Prompt as Administrator and run:
```
bcdedit /set hypervisorlaunchtype auto
```
Then restart.

**Issue: Antivirus blocks the installer**
Solution: Temporarily disable antivirus during installation only. Docker is safe software
from a reputable vendor.

---

## 4. Pulling the MOOSE Image

Once Docker Desktop is running, open a terminal (Git Bash recommended) and pull the
official MOOSE image:

```bash
docker pull idaholab/moose:latest
```

### What Happens During the Pull

You will see output like this:
```
latest: Pulling from idaholab/moose
6414378b6477: Pull complete
3d2430473443: Pull complete
...
Digest: sha256:a1b2c3d4...
Status: Downloaded newer image for idaholab/moose:latest
docker.io/idaholab/moose:latest
```

Each line represents one "layer" of the image being downloaded. Docker images are built in
layers — each layer is a set of filesystem changes (installing a package, copying files,
etc.). Layers are cached independently, so when Idaho National Laboratory releases an
updated image that only changes the top layer, you only need to download the changed layer.

**Download size:** The MOOSE image is approximately 2-4 GB. Download time depends on
your internet connection:
- 100 Mbps connection: approximately 5-10 minutes
- 25 Mbps connection: approximately 20-40 minutes

The download only needs to happen once. After that, the image is stored locally and
`docker run` starts instantly.

### What is Inside the Image

The `idaholab/moose:latest` image contains a complete, pre-built MOOSE environment:

| Component | Details |
|-----------|---------|
| Base OS | Ubuntu 22.04 LTS |
| Compiler | GCC with C++17 support |
| MPI | MPICH (optimized for MOOSE) |
| PETSc | Configured with MOOSE's required options |
| libMesh | Built against the MOOSE-compatible PETSc |
| HDF5 | For binary I/O |
| Python 3 | Including packages needed by TestHarness |
| MOOSE framework | Pre-compiled as `libmoose` |
| Test executable | `/opt/moose/bin/moose_test-opt` |
| Combined modules | `/opt/moose/bin/combined-opt` |

The two most important executables inside the container:

- **`/opt/moose/bin/moose_test-opt`** — The MOOSE test application. Includes all built-in
  kernels, boundary conditions, materials, and mesh generators. Suitable for running the
  Quick-Start examples and most tutorial problems. This is what you will use for the
  examples in this guide.

- **`/opt/moose/bin/combined-opt`** — The combined modules application. Includes all 25+
  optional physics modules (heat transfer, solid mechanics, Navier-Stokes, etc.) in
  addition to the framework. Use this for real physics problems that require specific
  module kernels.

### Available Image Tags

```bash
docker pull idaholab/moose:latest        # Latest stable release
docker pull idaholab/moose-dev:latest    # Development image (includes compilers)
```

The `moose-dev` image is larger but includes the full build toolchain, allowing you to
compile your own custom MOOSE applications inside the container.

### Verify the Image Downloaded Successfully

```bash
docker images
```

Expected output:
```
REPOSITORY         TAG       IMAGE ID       CREATED        SIZE
idaholab/moose     latest    a1b2c3d4e5f6   2 weeks ago    3.8GB
hello-world        latest    d2c94e258dcb   10 months ago  13.3kB
```

---

## 5. Running MOOSE Simulations

### 5.1 Basic Command Structure

Here is the fundamental command to run a MOOSE simulation from Windows:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/moose_test-opt -i your_input.i'
```

Every part of this command has a specific purpose. The following sections break down each
flag and explain the reasoning.

### 5.2 Every Flag Explained

#### `MSYS_NO_PATHCONV=1`

This environment variable prefix is specific to Git Bash on Windows.

**The problem:** Git Bash is built on top of a compatibility layer called MSYS2, which
tries to be helpful by automatically converting Unix-style paths in your commands to
Windows-style paths. For example, if you type `/bin/bash`, MSYS2 sees a Unix absolute
path and converts it to something like `C:/Program Files/Git/usr/bin/bash` — the
equivalent path on your Windows system.

This conversion is usually helpful when running native Windows programs from Git Bash. It
is disastrous when running Docker commands, because paths like `/bin/bash` and `/work`
are meant to refer to locations **inside the container** (which runs Linux), not on your
Windows machine. After MSYS2's "helpful" conversion, Docker receives a completely wrong
path that does not exist inside the container.

**The fix:** Setting `MSYS_NO_PATHCONV=1` tells MSYS2 to leave all paths exactly as you
typed them. Docker then receives the correct Linux paths for the container interior.

**Alternative approaches:**
- Use PowerShell instead of Git Bash (PowerShell does not mangle paths)
- Use double slashes: `//bin//bash` and `//work` (MSYS2 does not convert double-slash paths)
- Use Windows Command Prompt (cmd.exe) — also does not mangle paths

```bash
# Git Bash WITHOUT MSYS_NO_PATHCONV (broken — MSYS2 converts /bin/bash):
docker run --entrypoint /bin/bash idaholab/moose:latest
# Docker actually receives: --entrypoint C:/Program Files/Git/usr/bin/bash (WRONG)

# Git Bash WITH MSYS_NO_PATHCONV (correct):
MSYS_NO_PATHCONV=1 docker run --entrypoint /bin/bash idaholab/moose:latest
# Docker receives: --entrypoint /bin/bash (correct — refers to inside the container)
```

#### `docker run`

Creates a new container from the specified image and runs a command inside it. Every time
you run this, a fresh container is created from the clean image. No state is carried over
from previous runs (unless you explicitly use persistent volumes or `docker commit`).

#### `--rm`

Automatically removes the container after it exits.

Without `--rm`, stopped containers accumulate on your system. Each stopped container
takes up a small amount of disk space (just the filesystem diff from the image, typically
a few MB). Over many runs this adds up. `--rm` keeps your system clean automatically.

To check for accumulated stopped containers:
```bash
docker ps -a
```

To clean them up if you forgot `--rm`:
```bash
docker container prune
```

#### `-v "C:/Users/you/simulations:/work"`

This is the volume mount — the most critical flag for using Docker with MOOSE.

Format: `-v "host_path:container_path"`

- **`host_path`**: The folder on your Windows machine (use forward slashes, not backslashes)
- **`container_path`**: Where that folder appears inside the container

After mounting, any file in `C:/Users/you/simulations/` appears at `/work/` inside the
container, and vice versa. Files written to `/work/` by MOOSE appear in
`C:/Users/you/simulations/` on your Windows machine immediately.

**Path format rules for Git Bash:**
```bash
# Correct — use forward slashes and quote the whole -v value:
-v "C:/Users/you/simulations:/work"

# Also correct — use /c/ notation (Git Bash maps drives this way):
-v "/c/Users/you/simulations:/work"

# Wrong — backslashes will cause errors:
-v "C:\Users\you\simulations:/work"
```

#### `-w /work`

Sets the working directory inside the container. When the container starts, the shell's
current directory will be `/work/`. This means you can refer to input files by name
alone (e.g., `input.i`) rather than by full path (e.g., `/work/input.i`).

#### `--entrypoint /bin/bash`

Overrides the default command the container runs at startup. The `idaholab/moose` image
may have a default entrypoint configured. By setting it to `/bin/bash`, we tell Docker to
use the Bash shell as the entry point, allowing us to pass a shell command string via `-c`.

#### `idaholab/moose:latest`

The image to use. Format: `repository/name:tag`. The `latest` tag refers to the most
recently published version. You can pin to a specific version if needed (e.g.,
`idaholab/moose:2024.08`).

#### `-c '/opt/moose/bin/moose_test-opt -i your_input.i'`

The command string passed to bash. Bash executes this as a shell command inside the
container. The single quotes around the command are important — they prevent your Windows
shell from interpreting any special characters before passing the string to Docker.

### 5.3 Volume Mounting Explained

Volume mounting is the bridge between your Windows files and the container. Understanding
it clearly prevents most confusion when working with Docker and MOOSE.

```
  Windows filesystem (C:\)                  Docker container filesystem (/)
  ─────────────────────────                 ─────────────────────────────
  C:\
  └── Users\
      └── you\
          └── simulations\          ←──────── /work/
              ├── input.i           ←──mount──→ ├── input.i
              ├── sub_mesh.i        ←──mount──→ ├── sub_mesh.i
              ├── material.i        ←──mount──→ ├── material.i
              │                                  │
              │  (after simulation runs)         │  (MOOSE writes here)
              ├── input_out.e       ←──mount──→ ├── input_out.e
              ├── input_out.csv     ←──mount──→ ├── input_out.csv
              └── input_out.png     ←──mount──→ └── input_out.png
```

Key facts about volume mounts:
- The mount is **bidirectional** — changes flow both ways instantly
- The mount is **live** — files do not need to be copied in or out
- Files written inside the container to `/work/` appear on Windows immediately
- You can edit input files on Windows and run the simulation again without any additional
  steps
- The container does **not** have access to any other part of your Windows filesystem
  (only the mounted folder)

**Important:** Only mount the specific folder containing your simulation files. Do not
mount `C:/` (your entire C drive) as this gives the container unrestricted access to your
Windows filesystem.

### 5.4 Path Gotchas on Windows

This section collects all the path-related issues you are likely to encounter.

#### Issue 1: MSYS2 Path Conversion (Git Bash)

Already covered above. Solution: prefix commands with `MSYS_NO_PATHCONV=1`.

```bash
# Wrong (MSYS2 converts /bin/bash to a Windows path):
docker run --rm --entrypoint /bin/bash idaholab/moose:latest -c '...'

# Correct:
MSYS_NO_PATHCONV=1 docker run --rm --entrypoint /bin/bash idaholab/moose:latest -c '...'
```

#### Issue 2: Backslashes in Volume Mount Paths

Windows uses backslashes (`\`) in paths, but Docker expects forward slashes (`/`).

```bash
# Wrong — backslashes cause parsing errors:
-v "C:\Users\you\sims:/work"

# Correct — forward slashes only:
-v "C:/Users/you/sims:/work"

# Also correct — Git Bash drive notation:
-v "/c/Users/you/sims:/work"
```

#### Issue 3: Spaces in Windows Paths

If your path contains spaces, quote the entire `-v` value:

```bash
# Wrong — Docker cannot parse the path:
-v C:/Users/John Smith/simulations:/work

# Correct — quoted path with spaces:
-v "C:/Users/John Smith/simulations:/work"
```

#### Issue 4: Using PowerShell Instead of Git Bash

PowerShell does not mangle paths, so you do not need `MSYS_NO_PATHCONV`. However, the
line continuation character is a backtick (`` ` ``) instead of backslash (`\`):

```powershell
# PowerShell syntax:
docker run --rm `
  -v "C:/Users/you/simulations:/work" `
  -w /work `
  --entrypoint /bin/bash `
  idaholab/moose:latest `
  -c '/opt/moose/bin/moose_test-opt -i your_input.i'
```

#### Issue 5: WSL2 Paths

If you are running commands from inside a WSL2 session (not Git Bash), Windows drives are
available at `/mnt/c/`:

```bash
# WSL2 bash syntax — use /mnt/c/ instead of C:/:
docker run --rm \
  -v "/mnt/c/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/moose_test-opt -i your_input.i'
```

### 5.5 Running Examples

The following are complete, copy-paste commands. Replace `C:/Users/you/simulations` with
the actual path to the folder containing your input files.

#### Example 1: Running a Single Input File

Suppose you have `diffusion.i` in `C:/Users/you/simulations/`:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/moose_test-opt -i diffusion.i'
```

After the run, `diffusion_out.e` will appear in `C:/Users/you/simulations/`.

#### Example 2: Running a Parent/Sub-App Input (MultiApp)

MultiApp simulations involve a parent input file that launches sub-application instances.
Both files must be in the mounted volume. Suppose your files are:

```
C:/Users/you/simulations/
  parent.i
  sub.i
```

Run the parent only (it references `sub.i` internally):
```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/moose_test-opt -i parent.i'
```

The sub-application executable is also launched by MOOSE internally — you only invoke the
parent from the command line. Both `parent_out.e` and `sub_out.e` will appear in your
simulations folder after completion.

#### Example 3: Running with Parallel MPI

To run on multiple CPU cores simultaneously, use `mpirun`:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c 'mpirun -np 4 /opt/moose/bin/moose_test-opt -i diffusion.i'
```

The `-np 4` flag tells MPI to launch 4 parallel processes. MOOSE will automatically
partition the mesh across the 4 processes and run the solve in parallel. Results are
automatically recombined in the output file.

**Important:** Do not use more processes than the CPU cores you allocated to Docker
Desktop. Using `-np 8` when Docker has only 2 CPUs will actually slow things down.

**When is parallel useful?** For small tutorial problems, serial (no `mpirun`) is fine
and easier to debug. For large meshes (100k+ elements) or expensive physics, parallel
provides significant speedup.

#### Example 4: Interactive Container Session

Sometimes you want to explore the container, run multiple commands, or debug:

```bash
MSYS_NO_PATHCONV=1 docker run --rm -it \
  -v "C:/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest
```

The `-it` flags mean:
- `-i` (interactive) — keep stdin open so you can type commands
- `-t` (tty) — allocate a terminal so the output is formatted correctly

You will get a shell prompt inside the container:
```
root@a1b2c3d4e5f6:/work#
```

From here you can run multiple MOOSE commands, inspect the filesystem, check what
executables are available, and so on:

```bash
# Inside the container:
ls /opt/moose/bin/          # list available executables
/opt/moose/bin/moose_test-opt --help    # show MOOSE options
/opt/moose/bin/moose_test-opt -i diffusion.i    # run simulation
ls -lh                      # check output files
exit                        # leave the container
```

When you type `exit`, the container stops and is removed (because of `--rm`). Any files
you wrote to `/work/` remain on your Windows filesystem.

#### Example 5: Running All Quick-Start Cases in a Loop

If you have downloaded the 13 Quick-Start input files into a folder, you can run all of
them with a single command:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/you/quickstart:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c 'for i in case*.i; do
        echo "=== Running $i ===";
        /opt/moose/bin/moose_test-opt -i "$i" 2>&1;
        echo "=== Done with $i ===";
      done'
```

This iterates over every file matching `case*.i` in the mounted folder and runs each one.
Output files for each case will appear in `C:/Users/you/quickstart/` as they complete.

#### Example 6: Using the Combined Modules Application

For physics problems that require specific modules (heat transfer, solid mechanics, etc.),
use `combined-opt` instead of `moose_test-opt`:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i heat_conduction.i'
```

`combined-opt` includes everything in `moose_test-opt` plus all the physics module kernels
and materials. Use it whenever your input file references kernels from modules like
`HeatConduction`, `SolidMechanics`, `NavierStokes`, etc.

#### Example 7: Checking Available Kernels

Not sure which kernels are available in the container's executable? Use the `--show-type`
and `--registry` options:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt --registry 2>&1 | head -100'
```

This lists all registered object types (kernels, BCs, materials, etc.).

### 5.6 Checking Simulation Progress

MOOSE prints solve information to stdout, which Docker forwards to your terminal. A
typical solve looks like:

```
Framework Information:
MOOSE Version:          git commit ...
...

Mesh Information:
 Num Nodes:             121
 Num Elems:             100
...

Solving...
 0 Nonlinear |R| = 1.234e+00
      0 Linear |R| = 1.234e+00
      ...
      12 Linear |R| = 1.234e-08
 1 Nonlinear |R| = 2.345e-07
Converged!

Postprocessor Values:
...

Finished simulation in 0.23 seconds.
```

The solve is progressing normally if you see the "Nonlinear |R|" values decreasing each
iteration. If the values are increasing or the solve runs for many iterations without
converging, there may be a problem with the input file.

---

## 6. Docker for Custom MOOSE Applications

### 6.1 When You Need a Custom Application

The pre-built `moose_test-opt` and `combined-opt` executables cover a wide range of
physics problems. However, you need to build a custom application when you want to:

- Add your own physics (custom Kernels, Materials, BCs written in C++)
- Add custom postprocessors or outputs
- Couple to an external code
- Use MOOSE as the basis for your own specialized simulation tool

For custom applications, you need the `idaholab/moose-dev` image, which includes the
full compiler toolchain.

### 6.2 Pulling the Development Image

```bash
docker pull idaholab/moose-dev:latest
```

The dev image is larger (~5-8 GB) because it includes GCC, G++, `mpicxx`, `make`, CMake,
and all MOOSE header files and library archives needed for compilation.

### 6.3 Compiling a Custom Application Inside Docker

Suppose you have created a custom application using the `stork.sh` generator and the
source lives at `C:/Users/you/MyApp/`:

```
C:/Users/you/MyApp/
  Makefile
  src/
    base/
      MyAppApp.C
    kernels/
      MyKernel.C
  include/
    base/
      MyAppApp.h
    kernels/
      MyKernel.h
  test/
    ...
```

Start an interactive session with the dev image:

```bash
MSYS_NO_PATHCONV=1 docker run --rm -it \
  -v "C:/Users/you/MyApp:/myapp" \
  -w /myapp \
  --entrypoint /bin/bash \
  idaholab/moose-dev:latest
```

Inside the container, compile:

```bash
# Inside the container:
make -j4 METHOD=opt
```

The compiled executable (`myapp-opt`) appears in `/myapp/`, which means it also appears
in `C:/Users/you/MyApp/` on Windows. However, this executable is a Linux binary and can
only be run inside the container or on a real Linux system.

To run your compiled application on an input file:

```bash
# Still inside the container:
cd test
./run_tests -j4
# or run directly:
./myapp-opt -i path/to/input.i
```

### 6.4 Persistent Containers for Development

When actively developing a custom application, repeatedly typing the full `docker run`
command is tedious. Instead, start a named container and leave it running:

```bash
# Start a named, persistent container:
MSYS_NO_PATHCONV=1 docker run -d \
  --name myapp-dev \
  -v "C:/Users/you/MyApp:/myapp" \
  -v "C:/Users/you/simulations:/work" \
  -w /myapp \
  --entrypoint /bin/bash \
  idaholab/moose-dev:latest \
  -c 'while true; do sleep 60; done'
```

The `-d` flag runs the container in the background (detached). The `while true; do sleep 60; done` command keeps the container alive.

Now you can run commands inside the running container without creating a new one each time:

```bash
# Execute a command in the running container:
docker exec myapp-dev bash -c 'make -j4 METHOD=opt'

# Open an interactive shell in the running container:
docker exec -it myapp-dev bash
```

To stop and remove the container when you are done:

```bash
docker stop myapp-dev
docker rm myapp-dev
```

### 6.5 Saving a Compiled State with `docker commit`

If you want to save a container's state (including the compiled application) as a new
image that can be reused later:

```bash
# After compiling inside the container:
docker commit myapp-dev myapp:compiled
```

Now `myapp:compiled` is a new local image with your application already compiled. You can
run it later without recompiling:

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  myapp:compiled \
  -c '/myapp/myapp-opt -i input.i'
```

**Note:** `docker commit` creates a new local image. It does not push to Docker Hub. Use
`docker push` (after `docker tag` and creating a Docker Hub account) if you want to share
the image.

---

## 7. Troubleshooting

### 7.1 "Cannot connect to the Docker daemon"

**Full error:**
```
Cannot connect to the Docker daemon at unix:///var/run/docker.sock.
Is the docker daemon running?
```

**Cause:** Docker Desktop is not running.

**Fix:**
1. Check the system tray for the Docker whale icon. If it is not there, Docker Desktop
   is not running.
2. Open Docker Desktop from the Start menu.
3. Wait for the "Engine running" status (bottom left of the Docker Desktop window).
4. Try the command again.

If Docker Desktop is open but shows "Engine stopped":
```powershell
# PowerShell (as Administrator):
wsl --shutdown
```
Then click "Restart" in Docker Desktop.

### 7.2 "Error response from daemon: pull access denied" or "not found"

**Cause:** The image name is misspelled, or the image does not exist on Docker Hub.

**Fix:** Double-check the exact image name: `idaholab/moose:latest` (all lowercase). Try
pulling again:
```bash
docker pull idaholab/moose:latest
```

### 7.3 "The system cannot find the file specified" (Windows path error)

**Full error:**
```
docker: Error response from daemon: invalid volume specification:
'C:\Users\you\sims:/work': invalid mode: /work.
```

**Cause:** Backslashes in the volume path, or MSYS2 path conversion.

**Fix:**
```bash
# Use forward slashes:
-v "C:/Users/you/sims:/work"

# And prefix with MSYS_NO_PATHCONV=1 when using Git Bash:
MSYS_NO_PATHCONV=1 docker run --rm -v "C:/Users/you/sims:/work" ...
```

### 7.4 Output Files Not Appearing on Windows

**Symptom:** The simulation completes without errors, but you cannot find the output
`.e` or `.csv` files on Windows.

**Cause:** Either the volume mount path is wrong, or MOOSE is writing output to a
different directory.

**Diagnosis:**
```bash
# Run interactively and check where files end up:
MSYS_NO_PATHCONV=1 docker run --rm -it \
  -v "C:/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest

# Inside the container:
/opt/moose/bin/moose_test-opt -i input.i
ls -la    # check files in /work
```

If files appear in `/work/` inside the container, they should be on Windows. If not, the
`-v` mount is incorrect.

Check that the Windows path in the `-v` flag exactly matches the folder where your input
file is located. Use forward slashes and include the drive letter.

### 7.5 "Permission denied" on Mounted Files

**Symptom:** MOOSE cannot read your input file:
```
terminate called after throwing an instance of 'std::runtime_error'
what(): Unable to open file: input.i
```

**Cause:** Docker Desktop's file sharing is not configured to allow access to that folder.

**Fix for Hyper-V backend:**
1. Open Docker Desktop Settings
2. Resources → File Sharing
3. Add `C:\` (or the specific folder)
4. Apply & Restart

**Fix for WSL2 backend:** WSL2 should have full access to all Windows drives by default.
Restart Docker Desktop if the issue persists.

### 7.6 "No space left on device" Inside the Container

**Symptom:**
```
Error: No space left on device
```

**Cause:** Docker's virtual disk is full.

**Fix:**
```bash
# Remove unused containers, images, networks, and build cache:
docker system prune

# More aggressive — also remove unused images:
docker system prune -a

# Check how much space Docker is using:
docker system df
```

Note: `docker system prune` does not remove your running containers or named images that
are actively used.

### 7.7 Container Exits Immediately Without Output

**Symptom:** `docker run` returns with exit code 1 immediately, no output.

**Cause:** Usually a syntax error in the `-c` command string.

**Debug approach:**
```bash
# Remove the -c flag entirely and run interactively to see what happens:
MSYS_NO_PATHCONV=1 docker run --rm -it \
  -v "C:/Users/you/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest
```

Then type the command manually inside the container to see the exact error.

### 7.8 MOOSE Simulation Fails with "No registered object named X"

**Symptom:**
```
*** ERROR ***
No registered object named 'HeatConduction' in 'Kernels'.
```

**Cause:** You are using `moose_test-opt` but your input file requires a physics module
kernel that is only in `combined-opt`.

**Fix:** Switch to `combined-opt`:
```bash
-c '/opt/moose/bin/combined-opt -i your_input.i'
```

### 7.9 Very Slow Simulation Performance

**Possible causes and fixes:**

1. **Not enough memory allocated to Docker:**
   Docker Desktop Settings → Resources → Memory → Increase to 8+ GB → Apply & Restart.

2. **Not enough CPU cores allocated:**
   Docker Desktop Settings → Resources → CPUs → Increase → Apply & Restart.

3. **Volume mount performance (for disk-intensive simulations):**
   On Windows, Docker volume mounts between the Windows filesystem and Linux containers
   have some overhead. For very I/O-intensive simulations, consider placing files inside
   the Docker VM's native filesystem by using a named Docker volume instead:
   ```bash
   docker volume create moose-work
   docker run --rm -v moose-work:/work ...
   ```
   Copy files in with `docker cp`.

4. **Running too many MPI processes:**
   Use `-np` equal to the number of Docker CPUs, not more.

### 7.10 "Mounts denied: The path is not shared from the host"

**Full error:**
```
Error response from daemon: Mounts denied:
The path C:\Users\you\simulations
is not shared from OS X / Windows and is not known to Docker.
```

**Fix:**
- Docker Desktop → Settings → Resources → File Sharing
- Add the path that is being denied
- Apply & Restart Docker Desktop

### 7.11 Interactive Container Gets "^C" Stuck or Cannot Exit

If `Ctrl+C` does not stop a running simulation and you are stuck in an interactive
container, use `Ctrl+P` then `Ctrl+Q` to detach from the container without stopping it.
Then:

```bash
docker ps                   # find the container ID
docker stop <container_id>  # stop it
```

### 7.12 "WARNING: The requested image's platform does not match the host platform"

**Symptom:**
```
WARNING: The requested image's platform (linux/amd64) does not match
the host platform (linux/arm64/v8)
```

**Cause:** You are on an ARM-based machine (e.g., a newer laptop with Qualcomm Snapdragon
or Apple Silicon running in emulation). The MOOSE image is built for x86-64.

**Fix:** This is a warning, not an error. The container usually still works, but runs
under emulation and will be slower. Check the Docker Hub page for `idaholab/moose` to see
if an ARM image is available. If not, this is a known limitation.

---

## 8. Docker Commands Reference

### Essential Commands

| Command | Purpose | Example |
|---------|---------|---------|
| `docker pull <image>` | Download image from Docker Hub | `docker pull idaholab/moose:latest` |
| `docker run <image>` | Create and start a container | `docker run --rm idaholab/moose:latest` |
| `docker run -it <image>` | Start container with interactive terminal | `docker run -it idaholab/moose:latest bash` |
| `docker run -d <image>` | Start container in background | `docker run -d --name dev idaholab/moose-dev:latest` |
| `docker ps` | List running containers | `docker ps` |
| `docker ps -a` | List all containers (including stopped) | `docker ps -a` |
| `docker images` | List downloaded images | `docker images` |
| `docker exec <id> <cmd>` | Run command in existing container | `docker exec dev bash -c 'make'` |
| `docker stop <id>` | Stop a running container | `docker stop myapp-dev` |
| `docker rm <id>` | Remove a stopped container | `docker rm myapp-dev` |
| `docker rmi <image>` | Remove a downloaded image | `docker rmi idaholab/moose:latest` |
| `docker logs <id>` | View container output | `docker logs myapp-dev` |
| `docker cp <src> <dst>` | Copy files to/from container | `docker cp input.i dev:/work/` |
| `docker commit <id> <name>` | Save container state as image | `docker commit dev myapp:v1` |
| `docker system prune` | Remove unused containers and images | `docker system prune -a` |
| `docker system df` | Show disk usage by Docker | `docker system df` |
| `docker inspect <id>` | Show detailed container/image info | `docker inspect myapp-dev` |
| `docker volume create` | Create a named volume | `docker volume create moose-work` |
| `docker volume ls` | List named volumes | `docker volume ls` |
| `docker volume rm` | Remove a named volume | `docker volume rm moose-work` |

### Volume Mount Syntax Reference

```bash
# Basic mount — host:container:
-v "C:/path:/work"

# Read-only mount (container cannot modify files):
-v "C:/path:/work:ro"

# Named volume (Docker-managed, not a Windows path):
-v moose-data:/work

# Multiple mounts:
-v "C:/simulations:/work" -v "C:/configs:/config"
```

### Useful Flags for `docker run`

| Flag | Purpose |
|------|---------|
| `--rm` | Delete container on exit |
| `-it` | Interactive terminal |
| `-d` | Detached (background) |
| `-v host:container` | Mount volume |
| `-w /path` | Set working directory |
| `-e VAR=value` | Set environment variable |
| `--name myname` | Give container a name |
| `--cpus 4` | Limit CPU usage |
| `--memory 8g` | Limit memory usage |
| `--entrypoint /bin/bash` | Override default entrypoint |
| `-u $(id -u)` | Run as current user (Linux host only) |

### Checking Image Contents

```bash
# List files in a specific directory inside the image:
MSYS_NO_PATHCONV=1 docker run --rm \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c 'ls -la /opt/moose/bin/'

# Check MOOSE version inside the image:
MSYS_NO_PATHCONV=1 docker run --rm \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/moose_test-opt --version'

# Show all environment variables in the container:
MSYS_NO_PATHCONV=1 docker run --rm \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c 'env | sort'
```

---

## 9. Comparison: Docker vs WSL2 vs Native Linux

Choosing the right approach depends on your situation. Here is a detailed comparison:

| Feature | Docker | WSL2 | Native Linux |
|---------|--------|------|-------------|
| **Setup difficulty** | Easy — install Docker Desktop, pull image | Medium — enable WSL2, install Linux distro, compile MOOSE (~1 hour) | Hard — full OS install on machine or partition |
| **Time to first simulation** | 15-30 minutes (install + image pull) | 1-3 hours (install + MOOSE compilation) | Several hours (OS install + MOOSE compilation) |
| **Compilation required** | None (for standard simulations) | Yes — full MOOSE source build | Yes — full MOOSE source build |
| **Custom app development** | Possible (use moose-dev image) | Full support | Full support |
| **MPI performance** | Good (slight overhead from containerization) | Better (near-native Linux performance) | Best (no overhead) |
| **Disk I/O performance** | Good for container volumes; moderate for Windows mounts | Very good — WSL2 ext4 filesystem | Best |
| **Memory overhead** | Low (~200 MB for Docker engine) | Low (~200 MB for WSL2 VM) | None |
| **Disk usage** | ~4 GB per image (but images are layered and deduplicated) | ~10-20 GB for build environment and sources | ~20-30 GB |
| **File sharing with Windows** | Volume mounts (bidirectional, but path format matters) | Full access — `/mnt/c/` maps to `C:\` | Not applicable |
| **Windows IDE integration** | Good — edit files on Windows, run in container | Excellent — VS Code WSL extension runs natively in WSL2 | Not applicable |
| **Reproducibility** | Excellent — the same image runs identically everywhere | Good — tied to your specific build configuration | Depends on system configuration |
| **Sharing environment** | Trivial — share the image name | Harder — need to document build steps | Harder |
| **Internet required after setup** | No — image runs offline | No | No |
| **Version management** | Easy — pull different image tags | Manual — maintain separate builds | Manual |
| **HPC cluster use** | Singularity/Apptainer can convert Docker images | Cannot use directly on most HPC systems | Direct use |
| **Multiple MOOSE versions** | Easy — pull different tags, run side by side | Must manage separately compiled versions | Must manage separately compiled versions |
| **Recommended for** | Beginners, Windows users, reproducible workflows, sharing | Windows users who need full build capability | Active framework developers, production HPC runs |

### When to Choose Docker

Choose Docker if you:
- Are new to MOOSE and want to start running simulations today
- Are on Windows without 2+ hours to spend on compilation
- Need a reproducible, shareable simulation environment
- Are running standard physics from the built-in framework and combined modules
- Are teaching or running a course and need all students to have identical environments
- Want to try a new MOOSE version without affecting your existing installation

### When to Choose WSL2

Choose WSL2 if you:
- Need to develop custom C++ MOOSE objects and want a native-feeling Linux environment
- Work with VS Code and want the WSL2 extension for seamless editing
- Need better file I/O performance for disk-intensive simulations
- Plan to spend significant time with MOOSE and want a more permanent setup
- Are comfortable with Linux command-line environments

### When to Choose Native Linux

Choose Native Linux if you:
- Run MOOSE on an HPC cluster (native Linux is the only option)
- Are a framework developer contributing to MOOSE's C++ core
- Need maximum performance for production research simulations
- Already use Linux as your primary operating system

### HPC Compatibility Note

A significant advantage of the Docker approach is that Docker images can be converted to
Singularity/Apptainer container format, which is the containerization tool used on most
HPC clusters. If you develop and test simulations with the `idaholab/moose` Docker image
locally, you can often run the same simulation on a cluster by converting the image:

```bash
# On an HPC cluster with Singularity installed:
singularity pull moose.sif docker://idaholab/moose:latest
singularity exec --bind /scratch/you/project:/work moose.sif \
  /opt/moose/bin/combined-opt -i simulation.i
```

This makes Docker a particularly good choice if you will eventually move to cluster
computing.

---

## Quick Reference Card

Save this section for daily use.

### Start a MOOSE Simulation (Git Bash on Windows)

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/YOUR_USERNAME/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/moose_test-opt -i YOUR_INPUT.i'
```

### Start a MOOSE Simulation (PowerShell on Windows)

```powershell
docker run --rm `
  -v "C:/Users/YOUR_USERNAME/simulations:/work" `
  -w /work `
  --entrypoint /bin/bash `
  idaholab/moose:latest `
  -c '/opt/moose/bin/moose_test-opt -i YOUR_INPUT.i'
```

### Open an Interactive Shell

```bash
MSYS_NO_PATHCONV=1 docker run --rm -it \
  -v "C:/Users/YOUR_USERNAME/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest
```

### Run with 4 MPI Processes

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/YOUR_USERNAME/simulations:/work" \
  -w /work \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c 'mpirun -np 4 /opt/moose/bin/combined-opt -i YOUR_INPUT.i'
```

### Common Substitutions

| Replace this | With this |
|-------------|----------|
| `YOUR_USERNAME` | Your actual Windows username |
| `C:/Users/YOUR_USERNAME/simulations` | Path to your input file folder (forward slashes) |
| `YOUR_INPUT.i` | The name of your input file |
| `moose_test-opt` | `combined-opt` if your input uses physics module kernels |
| `-np 4` | Number of CPU cores to use (match Docker CPU allocation) |

---

*This guide covers Docker Desktop version 4.x and the `idaholab/moose:latest` image as
available in early 2026. Docker and MOOSE are actively developed; some details may change
in future releases. For the authoritative MOOSE documentation and module reference, see
the other documents in the `docs/` directory.*
