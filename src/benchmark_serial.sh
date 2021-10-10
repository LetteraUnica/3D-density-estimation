rm -rf density.bin
./generate_points.x 1000000
/usr/bin/time ./density_serial.x 80 0.127 points.bin