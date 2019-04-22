# [gls.py](gls.py) - Generalised Lomb-Scargle periodogram

gls.py is a standalone script for python2 and python3. It also runs from the command line.

Here is an example for a working plot:

<img src="demo_gls.png" width="50%"/>

### Installation

Create a folder and run there
```bash
wget https://raw.githubusercontent.com/mzechmeister/GLS/master/python/gls.py
python -c "import setuptools; setuptools.setup(name='gls')" develop --user
```


### Usage from shell command line
Run the demo with
```bash
python gls.py
```
or
```bash
python -m gls
```

A shortcut `gls` is very convenient
```bash
ln gls.py ~/bin/gls
```
Then you can use it for instance as
```bash
gls your_data_file -fend 1.
```

### Usage in Python
```python
from gls import Gls
```
and then see the example in
[gls.py#L104-L126](https://github.com/mzechmeister/GLS/blob/ea66a022e40208920bfeff150f6d461a748d4df3/python/gls.py#L104-L126)
