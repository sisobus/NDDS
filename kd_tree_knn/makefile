EXEC=./kdtree
CC       = g++
BOOST_INCLUDEDIR=/usr/local/boost_1_55_0
CCOPTS   = -c -g -I${BOOST_INCLUDEDIR} #-Wall
# -fpermissive
LINK     = g++ 
LINKOPTS = -pg -static -lm

.cpp.o:
	$(CC) $(CCOPTS) $<

all: kdtree

Node.o: Node.cpp Node.h QueryBox.h

QueryBox.o: QueryBox.cpp QueryBox.h

KDTree.o: Node.h KDTree.cpp KDTree.h QueryBox.h

UserInterface.o: UserInterface.h UserInterface.cpp Node.h KDTree.h QueryBox.h

UserInterfaceRandom.o:UserInterface.h UserInterfaceRandom.cpp Node.h KDTree.h QueryBox.h

kdtree: Node.o KDTree.o QueryBox.o UserInterface.o
	g++ -o ${EXEC} Node.o KDTree.o QueryBox.o UserInterface.o

kdtreerandom: Node.o KDTree.o QueryBox.o UserInterfaceRandom.o
	g++ -o ${EXEC}random Node.o KDTree.o QueryBox.o UserInterfaceRandom.o

clean:
	rm -f Node.o KDTree.o QueryBox.o UserInterface.o UserInterfaceRandom.o ${EXEC}

