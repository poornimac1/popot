#!/bin/bash

g++ -o main main.c `pkg-config --libs --cflags gsl` -O3