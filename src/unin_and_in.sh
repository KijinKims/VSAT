pip uninstall -y vsat
rm -rf vsat.egg-info
rm -rf dist
python setup.py sdist
pip install dist/vsat-1.0.tar.gz
