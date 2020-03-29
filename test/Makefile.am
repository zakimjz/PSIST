# Copyright (C) 2002 M.J. Zaki <zaki@badr.cs.rpi.edu> 
# Written by palmep2@cs.rpi.edu
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

INCLUDES = -I../src -I/usr/include/

# SPATIAL_INDEX_LIB_PATH = ../lib/spatialindex/lib
# MTREE_LIB_PATH = ../lib/M-tree
# GIST_LIB_PATH =  ../lib/M-tree/GiST/
# TOOLS_LIB_PATH = ../lib/tools/lib

CC = g++
AM_CXXFLAGS = -DUNIX -DPRINTING_OBJECTS

bin_PROGRAMS = query
#bin_PROGRAMS = bbindex query
#bin_PROGRAMS = gh_index gh_query

#bin_PROGRAMS = gh_index gh_query bbindex query

# bin_PROGRAMS = mtreeMain  mtreeQuery
# bin_PROGRAMS = tovector	

#tovector_SOURCES = toVector.cpp
#mtreeMain_SOURCES = mtreeMain.cpp
#mtreeQuery_SOURCES = mtreeQuery.cpp
#bbindex_SOURCES = bbindex.cpp
query_SOURCES = query.cpp

#gh_index_SOURCES = gh_index.cpp
#gh_query_SOURCES = gh_query.cpp


LDADD = -lpthread

# LDADD = -L$(GIST_LIB_PATH) -lgist -lrtree
# LDADD = ../src/libbbindex.a