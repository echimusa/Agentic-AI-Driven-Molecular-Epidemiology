# Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| `Permission denied` when running `./scripts/setup.sh` | The execute bit was lost (a zip file or a Windows checkout does this). Run `bash scripts/setup.sh ...` or `chmod +x scripts/*.sh`. |
| `/usr/bin/env: 'bash\r': No such file or directory` | The scripts were checked out with Windows line endings. Run `git config core.autocrlf false`, check the repository out again (`.gitattributes` sets LF for `.sh`), or `sed -i 's/\r$//' scripts/*.sh`. |
| `pwsh: command not found` | You are using the PowerShell scripts without PowerShell 7. Either use the Bash scripts (`scripts/*.sh`) or install `pwsh` (last section of [INSTALL-Ubuntu.md](INSTALL-Ubuntu.md)). |
| "running scripts is disabled on this system" (Windows) | `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`, then run the script again in the same window. |
| `setup.ps1` says Python or Rscript was not found right after `winget install` | Windows does not update the path of a running window. Open a new PowerShell window and rerun. |
| `error: externally-managed-environment` (pip on Ubuntu) | Do not install into the system Python. `setup.ps1` uses `.venv`; if you run pip yourself, use `.venv/bin/python -m pip`. |
| `pip install mappy` fails: `zlib.h` or a compiler is missing | `sudo apt-get install -y build-essential zlib1g-dev python3-dev`. On Windows `mappy` has no wheel: use WSL2 for Chapters 6 and 7. |
| An R package fails to compile | Install the system libraries named in the error (`sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev`) and run `Rscript r-packages.R <package>` again. |
| `check-install.ps1` marks a package DIFFERENT | The installed version is not the tested one. Results usually still agree to the printed precision; install the pinned version with `.venv/bin/python -m pip install name==version` if a `-Full` comparison fails. |
| `-Full` reports "reproduced byte for byte FAIL" for a table | A different version of IQ-TREE, TreeTime or a numerical library can change the last digits. Compare the two tables; the quick check (without `-Full`) is the test that the chapter's numbers hold. Chapter 5 is the most sensitive to the IQ-TREE version (tested with 2.0.7). |
| Chapter 5 fails with "iqtree2: not found" | `sudo apt-get install -y iqtree fasttree` (or `setup.ps1 -Chapters 5`). On Windows use WSL2. |
| `MemoryError` or the process is killed | The heaviest check (Chapter 9) peaked at 2.5 GB. Close other programs, or run chapters one at a time. |
| A script cannot download a file | The quick check never downloads. Rebuilding needs a network connection and, for NCBI, EBI and the UKHSA API, the sites to be reachable; the data change over time. |
| Garbled characters in the terminal (Windows) | The runner sets UTF-8 for Python. In an old console use Windows Terminal, or run `chcp 65001`. |
