ROOT_FLAGS = `root-config --cflags`
ROOT_LIBS  = `root-config --libs`
CFLAGS		 = -I.
RCXX = $(CFLAGS) $(ROOT_FLAGS) $(ROOT_LIBS) $(LIBS)

INC_PATH = -I/usr/include -Iinc
LIB_PATH = -L/usr/lib -L/usr/local/lib 
LIBS     = -lMinuit -lpthread
#LIBS     = -lboost_system -lMinuit -lpthread
           

# compiler flags:
#    -g    adds debugging information to the executable file
#    -Wall turns on most, but not all, compiler warnings

DBG    = -g -O3 -Wall 

CC     = g++ $(DBG) $(LIB_PATH) $(LIBS) $(INC_PATH) -lpthread -std=c++11
MIN_CC = g++ $(DBG) $(INC_PATH) $(CFLAGS) $(ROOT_FLAGS) -std=c++11

BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src
INC_DIR = inc

all: MonteCarlo

MonteCarlo: $(SRC_DIR)/MonteCarlo.cpp
	$(CC) -o $(BIN_DIR)/MonteCarlo $(SRC_DIR)/MonteCarlo.cpp $(RCXX)

clean:
	rm -f $(OBJ_DIR)/*.o $(BIN_DIR)/*
