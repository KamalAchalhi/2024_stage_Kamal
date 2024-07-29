import ants

tx = ants.new_ants_transform(dimension=2)

tx.set_parameters((0.9,0,0,1.1,10,11))

ants.write_transform(tx, '/envau/work/meca/users/2024_Kamal/output/output_script1/tx.mat')

tx2 = ants.read_transform('/envau/work/meca/users/2024_Kamal/output/output_script1/tx.mat')