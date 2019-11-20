import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():

    if (len(sys.argv) != 5):
        print("==========================================================")
        print("Usage:> %s <input_folder> <cable_length> <dx> <BCL>" % (sys.argv[0]))
        print("==========================================================")
        sys.exit(1)
    else:
        input_folder = sys.argv[1]
        cable_length = sys.argv[2]
        dx = sys.argv[3]
        bcl = sys.argv[4]

        offset = int((float(cable_length) / float(dx)) / 5)

        input_file_0 = input_folder + "/" + cable_length + "cm/" + bcl + "ms/sv-%d.dat" % (offset*0)
        input_file_1 = input_folder + "/" + cable_length + "cm/" + bcl + "ms/sv-%d.dat" % (offset*1)
        input_file_2 = input_folder + "/" + cable_length + "cm/" + bcl + "ms/sv-%d.dat" % (offset*2)
        input_file_3 = input_folder + "/" + cable_length + "cm/" + bcl + "ms/sv-%d.dat" % (offset*3)
        input_file_4 = input_folder + "/" + cable_length + "cm/" + bcl + "ms/sv-%d.dat" % (offset*4)

        data_0 = np.genfromtxt(input_file_0)
        data_1 = np.genfromtxt(input_file_1)
        data_2 = np.genfromtxt(input_file_2)
        data_3 = np.genfromtxt(input_file_3)
        data_4 = np.genfromtxt(input_file_4)

        # AP 1
        start = 0
        end = int(bcl)

        data_range_0 = data_0[start:end]
        vms_0 = data_range_0[:,1].tolist()
        max_vm_0 = max(vms_0)                             
        max_vm_index_0 = vms_0.index(max_vm_0) + start    
        start_0 = start
        end_0 = end

        data_range_1 = data_1[start:end]
        vms_1 = data_range_1[:,1].tolist()
        max_vm_1 = max(vms_1)
        max_vm_index_1 = vms_1.index(max_vm_1) + start
        start_1 = max_vm_index_1 - max_vm_index_0
        end_1 = start_1 + int(bcl)

        data_range_2 = data_2[start:end]
        vms_2 = data_range_2[:,1].tolist()
        max_vm_2 = max(vms_2)
        max_vm_index_2 = vms_2.index(max_vm_2) + start
        start_2 = max_vm_index_2 - max_vm_index_0
        end_2 = start_2 + int(bcl)

        data_range_3 = data_3[start:end]
        vms_3 = data_range_3[:,1].tolist()
        max_vm_3 = max(vms_3)
        max_vm_index_3 = vms_3.index(max_vm_3) + start
        start_3 = max_vm_index_3 - max_vm_index_0
        end_3 = start_3 + int(bcl)

        data_range_4 = data_4[start:end]
        vms_4 = data_range_4[:,1].tolist()
        max_vm_4 = max(vms_4)
        max_vm_index_4 = vms_4.index(max_vm_4) + start
        start_4 = max_vm_index_4 - max_vm_index_0
        end_4 = start_4 + int(bcl)

        dx = float(dx)
        delta_s1 = offset*dx
        delta_s2 = offset*dx
        delta_s3 = offset*dx

        delta_t1 = max_vm_index_2 - max_vm_index_1
        delta_t2 = max_vm_index_3 - max_vm_index_2
        delta_t3 = max_vm_index_4 - max_vm_index_3

        v1_ap_1 = float(delta_s1 / delta_t1)
        v2_ap_1 = float(delta_s2 / delta_t2)
        v3_ap_1 = float(delta_s3 / delta_t3)

        # Move to the second Action Potential
        start_0 = start_0 + int(bcl) + 1
        end_0 = start_0 + int(bcl)
        start_1 = start_1 + int(bcl) + 1
        end_1 = start_1 + int(bcl)
        start_2 = start_2 + int(bcl) + 1
        end_2 = start_2 + int(bcl)
        start_3 = start_3 + int(bcl) + 1
        end_3 = start_3 + int(bcl)
        start_4 = start_4 + int(bcl) + 1
        end_4 = start_4 + int(bcl)

        # AP 2
        data_range_0 = data_0[start_0:end_0]
        vms_0 = data_range_0[:,1].tolist()
        max_vm_0 = max(vms_0)                             
        max_vm_index_0 = vms_0.index(max_vm_0) + start_0

        data_range_1 = data_1[start_1:end_1]
        vms_1 = data_range_1[:,1].tolist()
        max_vm_1 = max(vms_1)                             
        max_vm_index_1 = vms_1.index(max_vm_1) + start_1

        data_range_2 = data_2[start_2:end_2]
        vms_2 = data_range_2[:,1].tolist()
        max_vm_2 = max(vms_2)                             
        max_vm_index_2 = vms_2.index(max_vm_2) + start_2

        data_range_3 = data_3[start_3:end_3]
        vms_3 = data_range_3[:,1].tolist()
        max_vm_3 = max(vms_3)                             
        max_vm_index_3 = vms_3.index(max_vm_3) + start_3

        data_range_4 = data_4[start_4:end_4]
        vms_4 = data_range_4[:,1].tolist()
        max_vm_4 = max(vms_4)                             
        max_vm_index_4 = vms_4.index(max_vm_4) + start_4

        delta_t1 = max_vm_index_2 - max_vm_index_1
        delta_t2 = max_vm_index_3 - max_vm_index_2
        delta_t3 = max_vm_index_4 - max_vm_index_3

        v1_ap_2 = float(delta_s1 / delta_t1)
        v2_ap_2 = float(delta_s2 / delta_t2)
        v3_ap_2 = float(delta_s3 / delta_t3)

        # Output the results
        output_filename = input_folder + "/" + cable_length + "cm/" + bcl + "ms/propagation-velocity.dat"
        output_file = open(output_filename,"w")

        output_file.write("%g %g %g %g\n" % (offset*dx,2.0*offset*dx,v1_ap_1*10.0,v1_ap_2*10.0))
        output_file.write("%g %g %g %g\n" % (2.0*offset*dx,3.0*offset*dx,v2_ap_1*10.0,v2_ap_2*10.0))
        output_file.write("%g %g %g %g\n" % (3.0*offset*dx,4.0*offset*dx,v3_ap_1*10.0,v3_ap_2*10.0))

        output_file.close()

        '''
        print("Cell 0")
        print("Max Vm = %g -- Index = %d" % (max_vm_0,max_vm_index_0))
        print("AP -- Start = %d -- End = %d" % (start_0,end_0))
        print("")

        print("Cell 100")
        print("Max Vm = %g -- Index = %d" % (max_vm_1,max_vm_index_1))
        print("AP -- Start = %d -- End = %d" % (start_1,end_1))
        print("")

        print("Cell 200")
        print("Max Vm = %g -- Index = %d" % (max_vm_2,max_vm_index_2))
        print("AP -- Start = %d -- End = %d" % (start_2,end_2))
        print("")

        print("Cell 300")
        print("Max Vm = %g -- Index = %d" % (max_vm_3,max_vm_index_3))
        print("AP -- Start = %d -- End = %d" % (start_3,end_3))
        print("")

        print("Cell 400")
        print("Max Vm = %g -- Index = %d" % (max_vm_4,max_vm_index_4))
        print("AP -- Start = %d -- End = %d" % (start_4,end_4))
        print("")
        '''
    

if __name__ == "__main__":
    main()
