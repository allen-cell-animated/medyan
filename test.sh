# Actin only
echo "actin only"
DIR=./actin_only
mkdir -p ./examples/${DIR}/out
./src/MEDYAN -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

# Linkers with actin
echo "high l over a"
DIR=./linker_actin/high_l_over_a
mkdir -p ./examples/${DIR}/out
./src/MEDYAN -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "low l over a"
DIR=./linker_actin/low_l_over_a
mkdir -p ./examples/${DIR}/out
./src/MEDYAN -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "similar l over a"
DIR=./linker_actin/similar_l_over_a
mkdir -p ./examples/${DIR}/out
./src/MEDYAN -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

# Motors with actin
echo "high m over a"
DIR=./motor_actin/high_m_over_a
mkdir -p ./examples/${DIR}/out
./src/MEDYAN -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "low m over a"
DIR=./motor_actin/low_m_over_a
mkdir -p ./examples/${DIR}/out
./src/MEDYAN -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "similar m over a"
DIR=./motor_actin/similar_m_over_a
mkdir -p ./examples/${DIR}/out
./src/MEDYAN -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

echo "extremly high m over a"
DIR=./motor_actin/extremly_high_m_over_a
mkdir -p ./examples/${DIR}/out
./src/MEDYAN -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

# Motor linker actin
echo "motor linker actin"
DIR=./motor_linker_actin
mkdir -p ./examples/${DIR}/out
./src/MEDYAN -s ./examples/${DIR}/systeminput.txt -i ./examples/${DIR} -o ./examples/${DIR}/out

