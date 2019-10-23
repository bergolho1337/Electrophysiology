import os
import sys
import numpy as np

def main():

	if len(sys.argv) != 1:
                print("-------------------------------------------------------------------------")
                print("Usage:> python %s " % sys.argv[0])
                print("-------------------------------------------------------------------------")
                return 1

        filepath = 'inputs/param_SS_1000_configs.txt'
        with open(filepath) as fp:
                #line = fp.readline()
                cnt = 0
                #while line:
                for line in fp:
                        #print("Line {}: {}".format(cnt, line.strip()))
                        #line = fp.readline()
                        if (len(line) > 0):
                                data = line.split('|')
                                params = data[1].split()
                                gates = data[4].split()

                                atpi = float(params[0])
                                Ko = float(params[1])
                                Ki = float(params[2])
                                Vm_modifier = float(params[3])
                                GNa_modifier = float(params[4])
                                GCaL_modifier = float(params[5])
                                INaCa_modifier = float(params[6])

                                sucess = int(data[2].split()[0])

                                sv = []
                                sv.append(float(data[3].split()[0]))

                                for i in range(11):
                                    sv.append(float(gates[i]))

                                if (sucess):
                                    output_filepath = "outputs/atpi_%g-Ko_%g-Ki_%g-Vm:mod_%g-GNa:mod_%g-GCaL:mod_%g-INaCa:mod_%g" % (atpi,Ko,Ki,Vm_modifier,GNa_modifier,GCaL_modifier,INaCa_modifier)
                                    output_file = open(output_filepath,'w')

                                    output_file.write("%g\n" % (atpi))
                                    output_file.write("%g\n" % (Ko))
                                    output_file.write("%g\n" % (Ki))
                                    output_file.write("%g\n" % (Vm_modifier))
                                    output_file.write("%g\n" % (GNa_modifier))
                                    output_file.write("%g\n" % (GCaL_modifier))
                                    output_file.write("%g\n" % (INaCa_modifier))

                                    for i in range(12):
                                        output_file.write("%g\n" % (sv[i]))

                                    output_file.close()

                                #print("atpi=%g | Ko=%g | Ki=%g | Vm_mod=%g | GNa_mod=%g | GCaL_mod=%g | INaCa_mod=%g" % (atpi,Ko,Ki,Vm_modifier,GNa_modifier,GCaL_modifier,INaCa_modifier))
                                #print("propagation sucess = %d" % sucess)
                                #for i in range(12):
                                #    print("sv[%d]=%g" % (i,sv[i]))
                                #print("\n")

                                cnt += 1

if __name__ == "__main__":
        main()
