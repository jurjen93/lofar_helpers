# Step by step for running a model
## Install requirements
```
python -m venv venv
source activate venv/bin/activate
pip install -r requirements.txt
```

## Convert fits files to npz (numpy compressed)
```python
python pre_processing_for_ml.py <path-to-fits-dir>
```

## (Optional) Copy files to /dev/shm for fast dataloading
```python
find <path-to-fits-dir> -type f -name "*.npz" | xargs -n 1 -P 8 -i rsync -R {} /dev/shm
```

## Run neural network training
```python
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
