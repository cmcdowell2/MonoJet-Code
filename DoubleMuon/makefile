rflags = -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs` -lFoam -lMathMore

analyze: MonoJets.C MonoJets.h MonoJetsInit.C
	g++ -Wno-deprecated MonoJets.C MonoJetsInit.C -o analyze $(rflags)
