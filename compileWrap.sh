cd src/main/java/edu/stevens/dhutchis/accumuloiter/
nvcc -Xcompiler -fPIC -o Wrap.so -shared -I$JAVA_HOME/include -I$JAVA_HOME/include/linux  *.cpp *.cu
cd -
