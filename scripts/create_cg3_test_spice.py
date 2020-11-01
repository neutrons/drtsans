# From /HFIR/CG3/IPTS-17241/exp549/Datafiles$ gvim BioSANS_exp549_scan0020_0001.xml
# To  BioSANS_exp549_scan0020_0001_test.xml
import os

# 1. set detector_trans pos="0.000":  Detecotr 0 and (192 * 256 - 1) are at symmetric position

spice_name = '/HFIR/CG3/IPTS-17241/exp549/Datafiles/BioSANS_exp549_scan0020_0001_test.xml'
with open(spice_name, 'r') as spice:
    spice_content = spice.read()

print(spice_content)
main_det_start_flag = '<Detector type="INT32[192,256]">'
main_det_end_flag = '</Detector>'

# part 1: beginning to <Detector type="INT32[192,256]">
spice_part_1 = spice_content.split(main_det_start_flag)[0]
# part 2: counts
main_det_count_str = spice_content.split(main_det_start_flag)[1].split(main_det_end_flag)[0]
# part 3: from '</Detector>' to the end
spice_part_rest = spice_content.split(main_det_end_flag, 1)[1]

# Now work on 'part 3'
wing_det_start_flag = '<DetectorWing type="INT32[160,256]">'
wing_det_end_flag = '</DetectorWing>'

spice_part_3 = spice_part_rest.split(wing_det_start_flag)[0]
wing_det_counts_str = spice_part_rest.split(wing_det_start_flag)[1].split(wing_det_end_flag)[0]
spice_part_4 = spice_part_rest.split(wing_det_start_flag)[1].split(wing_det_end_flag)[1]

# Set main detector counts
counts_array = main_det_count_str.split()
print(len(counts_array))
num_pixels = len(counts_array)

num_tags = 0
for iws in range(num_pixels):
    if (iws - 127) % 256 == 0:
        num_tags += 1
        tube_index = (iws - 127) // 256
        counts_array[iws] = str(tube_index + 1)

# form the count string again
new_counts_str = ''
for count in counts_array:
    new_counts_str += f'{count}  '

# set wing detectors counts
wing_count_array = wing_det_counts_str.split()
print(f'wing detector size: {len(wing_count_array)}')
num_pixels = len(wing_count_array)

num_tags = 0
for iws in range(num_pixels):
    if (iws - 127) % 256 == 0:
        num_tags += 1
        tube_index = (iws - 127) // 256
        wing_count_array[iws] = str((tube_index + 1) * 1000)

new_wing_counts_str = ''
for count in wing_count_array:
    new_wing_counts_str += f'{count}  '

# Create new SPICE file
new_spice = f'{spice_part_1}\n' \
            f'{main_det_start_flag}\n{new_counts_str}\n{main_det_end_flag}\n' \
            f'{spice_part_3}\n' \
            f'{wing_det_start_flag}\n{new_wing_counts_str}\n{wing_det_end_flag}\n' \
            f'{spice_part_4}'

# output
with open('new.xml', 'w') as out:
    out.write(new_spice)
print(os.getcwd())
