all: cpe.pdf

concentration-vs-distance.pdf impedance-vs-frequency.pdf potential-energy-vs-distance.pdf resistivity-vs-distance.pdf : RawDataGraphene.txt compute.py impedance.py
	python impedance.py

iterate.dat optimal-concentration.pdf optimal-potential.pdf optimal-resistivity-and-concentration.pdf rho.csv : RawDataGraphene.txt compute.py transposed_resistivity.py
	python transposed_resistivity.py

poisson.pdf resistivity.pdf sigma-vs-V.pdf : poisson.py
	python poisson.py

cpe.aux cpe.log cpe.pdf : concentration-vs-distance.pdf cpe.tex impedance-vs-frequency.pdf optimal-Z-fit.pdf optimal-potential.pdf optimal-resistivity-and-concentration.pdf poisson.pdf potential-energy-vs-distance.pdf resistivity-vs-distance.pdf resistivity.pdf sigma-vs-V.pdf
	pdflatex cpe.tex && pdflatex cpe.tex

