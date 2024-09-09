# Step by step for running a model
## Install requirements
```shell
python -m venv venv
source activate venv/bin/activate
pip install -r requirements.txt
```

## Get training data
```shell
mkdir <path-to-fits-dir>
cd <path-to-fits-dir>

ONLINE_DATA_PATH=https://public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/
wget -r --spider --no-parent $ONLINE_DATA_PATH 2>&1 | grep -o 'https://[^ ]*' | grep -E '^.*\/{1}[^\/]+?\.[^\/]+?$' > urls.txt
xargs -n 1 -P 16 wget -q -r -np --no-clobber < urls.txt 
```

## Preprocess fits files and convert to npz (numpy compressed)
This will put the converted `*.npz` files in the same dir as the input `*.fits` files.
```shell
python pre_processing_for_ml.py <path-to-fits-dir>
```

## (Optional) Copy files to /dev/shm for fast dataloading
```shell
find <path-to-fits-dir> -type f -name "*.npz" | xargs -n 1 -P 8 -i rsync -R {} /dev/shm
```

## Run neural network training
```shell
python train_nn.py /dev/shm/<fullpath>
```
The dataloader expect filetree to be in the following format:
```text
<filepath>
  |- continue
  |- stop
  |- continue_val
  |- stop_val
```

### Inference guide

Install this package as follows:
```
pip install -e "git+https://github.com/jurjen93/lofar_helpers.git#egg=cortex&subdirectory=neural_networks"
```

You can then initialize the predictor as follows:

```Python
from cortex.predictors import StopPredictor
predictor = StopPredictor(cache="/your/model/cache", device="cpu", model="name_of_the_model")
```

The model cache will be used to store the to-be-downloaded models. This will then reuse downloaded models on subsequent runs of the program.

```Python
predictor.predict(input_path="/your/input/file.fits")
```