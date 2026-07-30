# Tutorial data staging area

This folder is for data that are useful for tutorials but should not be committed to GitHub.

The intended public mirror is:

```text
https://popgen.dk/albrecht/open/tutorial_data/
```

Sync this folder to the public mirror from the repository root with:

```bash
rsync -a -P tutorial_data/ kelly.popgen.dk:/kellyData/home/albrecht/public/open/tutorial_data
```

Recommended layout:

```text
tutorial_data/
  admixture/
  ngsadmix/
  evaladmix/
```

Each subfolder should contain the exact files expected by the corresponding tutorial. The tutorial `README.md` should use local paths such as:

```bash
DATA_ROOT=${TUTORIAL_DIR}/../tutorial_data
DATA_DIR=${DATA_ROOT}/ngsadmix
```

and should also document the public URL:

```bash
DATA_URL=https://popgen.dk/albrecht/open/tutorial_data/ngsadmix
```

Large files in this directory are ignored by git. Keep only lightweight `README.md` files tracked.
