# Publishing the repository on GitHub

The repository https://github.com/echimusa/Agentic-AI-Driven-Molecular-Epidemiology is empty, with `main` as its default branch. To make the first push from the unzipped folder:

```bash
unzip Agentic-AI-Driven-Molecular-Epidemiology.zip
cd Agentic-AI-Driven-Molecular-Epidemiology
git init -b main
git add .
git commit -m "Code, data and checks for the 13 chapters"
git remote add origin https://github.com/echimusa/Agentic-AI-Driven-Molecular-Epidemiology.git
git config http.postBuffer 524288000      # the repository is about 150 MB; this avoids "RPC failed" over HTTPS
git push -u origin main
```

GitHub asks you to sign in: use a personal access token as the password over HTTPS, or `gh auth login`, or an SSH remote (`git@github.com:echimusa/Agentic-AI-Driven-Molecular-Epidemiology.git`).

The repository holds only code, data, figures and check scripts. The lecture decks and the Word text of the chapters are not included, and `.gitignore` keeps `*.pptx`, `*.docx`, `chapters/*/lecture/` and `chapters/*/manuscript/` out of any commit, so they cannot be added by accident if you copy them into the folder to run the text checks (`--docx`).

Three things to decide before you push:

1. **CARD / RGI output files** (`chapters/chapter09/data/rgi_outputs/`, 12 files, 5 MB). They contain CARD material, and CARD's terms allow free non-commercial use only by individuals at academic institutions; everyone else needs written permission from McMaster University (see [SOURCE-LICENCES.csv](SOURCE-LICENCES.csv) and [DATA.md](DATA.md)). To keep them out of the public repository, run this before `git add .`: `echo "chapters/chapter09/data/rgi_outputs/" >> .gitignore`. Chapter 9's two-tool comparison then needs those files from your own copy; the result tables stay in the repository.
2. **Data terms**: [DATA.md](DATA.md) records the terms as understood and is not legal advice.
3. **Licence files**: `LICENSE` (MIT, code) and `LICENSE-CC-BY-4.0.txt` (text and figures) are in place; see [LICENSING.md](../LICENSING.md).

No file is above 21 MB (the Chapter 10 genotype extract is stored compressed), so GitHub raises no size warning. The scripts keep their execute permission, so `./scripts/setup.sh` works after `git clone` on Ubuntu.
