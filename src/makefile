TOPTARGETS := all clean

# SUBDIRS := $(wildcard */.)
SUBDIRS := auto_simplex ga param_simplex

PACKDIR = ../package

$(TOPTARGETS): $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: $(TOPTARGETS) $(SUBDIRS)

package: 
#	rm -rf $(PACKDIR)
#	cp -a ../package_template $(PACKDIR)
	mkdir -p $(PACKDIR)/src/auto_simplex
	mkdir -p $(PACKDIR)/src/ga
	mkdir -p $(PACKDIR)/src/param_simplex
	cp -a exec_single.sh recipe.txt settings.inp flow.py peakSeparationForXPS.py make_result_figures_ppt2.py $(PACKDIR)/src
	cp -a ./auto_simplex/*.exe ./auto_simplex/*.py $(PACKDIR)/src/auto_simplex
	cp -a ./ga/*.exe ./ga/*.py  $(PACKDIR)/src/ga
	cp -a ./param_simplex/*.exe ./param_simplex/*.py $(PACKDIR)/src/param_simplex



#	cp -a $(TARGET) automatic_xps_peak_separation.py automatic_xps_peak_separation_single.py exec_all.sh make_result_figures_ppt2.py myfunc.py plot_result_shirley_and_linear.py recipe.txt settings_Linear.inp settings_Shirley.inp $(PACKDIR)

