# search Fortran compiler
FC = `whereis gfortran f77 g77 | sed '/:$$/d; s/:.*//; s/gfortran/& -fbackslash -w/; q' `
gp = `gnuplot --version`
#
GLS:
	@ echo $(FC) GLS.f -o GLS.e
	$(FC) GLS.f -o GLS.e
#	gfortran  GLS.f -o GLS.e
	@ (echo GLSDir=`pwd`; cat GLSbin)  >  ~/bin/GLS
	chmod u+x ~/bin/GLS
	@ echo $(gp)
	@ if [[ "$(gp)" >  "gnuplot 4.3" ]]; then\
	   sed -i "s:pause mouse key:pause mouse button2:" ~/bin/GLS;\
	elif [[ "${gp}" =~ "gnuplot 4.2" ]]; then\
	   echo "pause patch for Gnuplot version 4.2";\
	   sed -i "/^zoomfile=/ s,=.*,='"`pwd`"/zoom.gnu'," pausepatch.gnu; \
	   sed -i "s,gnuplot &,gnuplot -ctrlq \&," ~/bin/GLS;\
	   sed -i "s:pause mouse key:\"load '\$$GLSDir/pausepatch.gnu'\":" ~/bin/GLS;\
	fi

