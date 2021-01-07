## all         : default task
all:
	pip install -e .

test:
	python -m unittest -v tests/*.py