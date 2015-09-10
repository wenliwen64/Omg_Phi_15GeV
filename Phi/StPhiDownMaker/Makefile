# =============================== 
# Makefile for StrAnalyMaker 
# =============================== 
ROOTGLIBS     = `root-config --glibs`
ROOTCFLAGS    = `root-config --cflags` 

CXX           = g++ 
CXXFLAGS      = -g -Wall -fPIC
#CXXFLAGS      = -g -Wall -fPIC -pthread -m64 -I/home/lwen/Downloads/root/include# ->home
#CXXFLAGS      = -I/usr/local/include/root ->office
SOFLAGS       = -shared 

CXXFLAGS     += $(ROOTCFLAGS)

MYCLASSES = StPhiDownMaker.o StPhiDownMakerDict.o 

# =============================== 
StrAnalyMaker: StPhiDownMaker.cc 
# ------------------------------- 
	echo $(ROOTGLIBS)
	echo $(ROOTCFLAGS)
	$(CXX) $(CXXFLAGS) -c StPhiDownMaker.cc -o StPhiDownMaker.o 
	rootcint -f StPhiDownMakerDict.cc -c StPhiDownMaker.hh Linkdef.h
	$(CXX) $(CXXFLAGS) -c StPhiDownMakerDict.cc -o StPhiDownMakerDict.o 	
	$(CXX) $(SOFLAGS) $(MYCLASSES) $(ROOTGLIBS) -o StPhiDownMaker.so  
	#cp StPhiDownMaker.so /usr/local/lib/root

# =============================== 
install:
	cp StPhiDownMaker.so /usr/local/lib/root
clean: 
	rm -f StPhiDownMaker.o StPhiDownMakerDict.o StPhiDownMakerDict.cc StPhiDownMakerDict.h StPhiDownMaker.so StPhiDownMaker_cc.so StPhiDownMaker_cc.d 
# =============================== 
