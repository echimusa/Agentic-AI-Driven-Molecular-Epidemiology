# Install on Windows, using PowerShell

There are two routes. **Route A (WSL2) is the supported one and covers all 13 chapters.** Route B runs most chapters natively. **Neither route has been tested on a Windows machine**: the steps follow the documentation of `winget`, WSL2 and the packages, and the Ubuntu steps they rely on were tested (see the root README). If a step fails, [TROUBLESHOOTING.md](TROUBLESHOOTING.md) has the likely causes.

## Route A: WSL2 with Ubuntu (all chapters)

Most command-line tools of the book (IQ-TREE 2, MAFFT, FastTree, bcftools, PLINK, BEAST 2, AMRFinderPlus) and the `mappy` package have no ready Windows build, so the easy way is to run Ubuntu inside Windows.

1. Open PowerShell **as administrator** and run:

   ```powershell
   wsl --install -d Ubuntu-24.04
   ```

   Restart when asked, then open "Ubuntu 24.04" from the Start menu and create a Linux user.
2. Inside Ubuntu, follow [INSTALL-Ubuntu.md](INSTALL-Ubuntu.md) from step 1. It uses the Bash scripts (`scripts/setup.sh`, `scripts/run-chapter.sh`); PowerShell is not needed inside WSL2. Keep the repository in your Linux home folder (for example `~/Agentic-AI-Driven-Molecular-Epidemiology`), not under `/mnt/c`, which is much slower.
3. To open the results in Windows, use `\\wsl$\Ubuntu-24.04\home\<user>\Agentic-AI-Driven-Molecular-Epidemiology` in File Explorer.

WSL2 needs Windows 10 version 2004 or later, or Windows 11, with virtualisation enabled in the BIOS.

## Route B: native Windows (PowerShell 5.1 or 7)

Expected to work for Chapters 1 to 4 and 8 to 13. Chapter 5 needs IQ-TREE 2 and FastTree, and Chapters 6 and 7 need `mappy`, which has no Windows wheel (it compiles with the Visual C++ Build Tools): use Route A for those three.

1. Allow scripts for this session, and install Git:

   ```powershell
   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
   winget install --id Git.Git -e
   ```

   Open a **new** PowerShell window so that `git` is on the path.
2. Get the repository:

   ```powershell
   git clone https://github.com/echimusa/Agentic-AI-Driven-Molecular-Epidemiology.git
   cd Agentic-AI-Driven-Molecular-Epidemiology
   ```
3. Install Python 3.12, R, pandoc, the Python environment and the R packages:

   ```powershell
   .\scripts\setup.ps1 -Chapters 1,2,3,4,8,9,10,11,12,13
   ```

   `winget` installs Python 3.12, R and Pandoc. Windows does not refresh the path of a running window, so if the script says that Python or Rscript was not found, open a new PowerShell window and run it again. R packages come from CRAN as binaries and go to `R\library` inside the repository.
4. Check and run:

   ```powershell
   .\scripts\check-install.ps1 -Chapters 12,13
   .\scripts\run-chapter.ps1 -Chapter 13
   .\scripts\run-chapter.ps1 -Chapter 13 -Script ex13_5_sir_gi.py
   ```

Notes for Route B:

- The scripts read and write UTF-8; the runner sets `PYTHONUTF8=1`, so no change to the console is needed.
- If a path becomes too long, enable long paths: `git config --system core.longpaths true` (administrator) and turn on "Enable Win32 long paths" in Windows.
- Windows PowerShell 5.1 and PowerShell 7 both run the scripts; examples use `.\scripts\...`, which works in both.
- The reference results were produced on Ubuntu 24.04. A different operating system or library version can change the last digits of a result table, so a byte-for-byte difference in `-Full` is a prompt to compare the numbers, not proof of an error.
