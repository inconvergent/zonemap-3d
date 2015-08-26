
PYX = $(shell find ./src -iname "*.pyx"| sort)

.PHONY: pyx clean

all: clean pyx

pyx:
	python setup.py build_ext --inplace
	cython -a $(PYX)
	python setup.py install --user

clean:
	rm -rf *.pyc
	rm -rf build
	rm -rf dist
	rm -rf *.egg-info

