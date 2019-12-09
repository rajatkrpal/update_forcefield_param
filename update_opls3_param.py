#!/usr/bin/python

import os, shutil
import sys
import re
import sqlite3 as lite
import math


#parse atoms with its properties from the parameter file
def extract_atomtable(paramfile):
    Atom =[]
    row2 = True
    while row2:
        next_line = paramfile.next()
        row2 = next_line.split()
        if not row2:
            break
        else:
            Atom.append(row2)

    atomtable_column = zip(*Atom)
    number = atomtable_column[0]
    charge = atomtable_column[1]
    sigma = atomtable_column[2]
    epsilon = atomtable_column[3]
    #print Atom
    total_number_of_atoms = len(number)
    return number, charge, sigma, epsilon

#parse bond stretch parameters from the parameter file
def extract_stretch(paramfile):
    Stretch =[]
    row2 = True
    while row2:
        next_line = paramfile.next()
        row2 = next_line.split()
        if not row2:
            break
        else:
            Stretch.append(row2)
    
    stretch_column = zip(*Stretch)
    atom1_stretch = stretch_column[0]
    atom2_stretch = stretch_column[1]
    k_st = stretch_column[2]
    r_eq = stretch_column[3]
    return (atom1_stretch, atom2_stretch, k_st, r_eq)

#parse bond angles between atoms from the parameter file
def extract_bend(paramfile):
    Bend =[]
    row2 = True
    while row2:
        next_line = paramfile.next()
        row2 = next_line.split()
        if not row2:
            break
        else:
            Bend.append(row2)

    bend_column = zip(*Bend)
    bend_atom1 = bend_column[0]
    bend_atom2 = bend_column[1]
    bend_atom3 = bend_column[2]
    bend_k = bend_column[3]
    bend_theta_eq = bend_column[4]
    
    return (bend_atom1, bend_atom2, bend_atom3, bend_k, bend_theta_eq)

#parse the torsion parameters from the parameter file
def extract_torsion(paramfile):
    Torsion =[]
    row2 = True
    while row2:
        next_line = paramfile.next()
        row2 = next_line.split()
        if not row2:
            break
        else:
            Torsion.append(row2)

    torsion_column = zip(*Torsion)
    torsion_atom1 = torsion_column[0]
    torsion_atom2 = torsion_column[1]
    torsion_atom3 = torsion_column[2]
    torsion_atom4 = torsion_column[3]
    torsion1 = torsion_column[4]
    torsion2 = torsion_column[5]
    torsion3 = torsion_column[6]
    torsion4 = torsion_column[7]

    
    #print torsion_column

    return (torsion_atom1, torsion_atom2, torsion_atom3, torsion_atom4, torsion1, torsion2, torsion3, torsion4)

#parse the improper torsion parameters from the parameter file
def extract_improper_torsion(paramfile):
    Improper = []
    row2 = True
    while row2:
        next_line = paramfile.next()
        row2 = next_line.split()
        if not row2:
            break
        else:
            Improper.append(row2)
    
    improper_column = zip(*Improper)
    improper_atom1 = improper_column[0]
    improper_atom2 = improper_column[1]
    improper_atom3 = improper_column[2]
    improper_atom4 = improper_column[3]
    improper_torsion2 = improper_column[4]
    
    #print improper_atom4

    
    return (improper_atom1, improper_atom2, improper_atom3, improper_atom4, improper_torsion2)
   
def extract_nb14_scaling_factor(number, paramfile):
    total_number_of_atoms = len(number)
    zeros = [0 for n in range(0, total_number_of_atoms)]
    Nb14scale_coul=[{} for i in range(0, total_number_of_atoms+1)]
    Nb14scale_lj=[{} for j in range(0, total_number_of_atoms+1)]
    row2 = True
    while row2:
        next_line = paramfile.next()
        row2 = next_line.split()
        #print row2
        if not row2:
            break
        else:
            nb14atom1 = int(row2[0])
            nb14atom4 = int(row2[1])
            nb14coulomb = float(row2[2])
            nb14lj = float(row2[3])
            
            #print nb14atom1, nb14atom4
            Nb14scale_lj[nb14atom1][nb14atom4] = nb14lj
            Nb14scale_coul[nb14atom1][nb14atom4] = nb14coulomb
            
    return (Nb14scale_coul, Nb14scale_lj)
    
# adds atom information, charge, and nonbonded paramters to the particle table, update sigma and epsilon values in the nonbonded_param table.
#Also returns the atomnumber, charge, sigma, epsilon for updating the nonbonded 1-4 interaction tables.
def add_atomtable_to_dms(con, num, charge, sigma, epsilon): 
    #num, charge, sigma, epsilon = extract_atomtable(paramfile)
    with con:
        cur = con.cursor()
        del_nonbonded_param = "DELETE FROM nonbonded_param"
        cur.execute(del_nonbonded_param)
        for i in range(0,len(num)):
            index = int(num[i])
            chrg = float(charge[i])
            sig = float(sigma[i])
            epsi = float(epsilon[i])
            #print index , chrg
            
            
            nonbonded_param_cmd = """INSERT INTO nonbonded_param(sigma, epsilon)  
            SELECT {0} , {1} 
            EXCEPT 
            SELECT sigma, epsilon from nonbonded_param WHERE sigma = {0} AND epsilon = {1}""" .format(sig, epsi)
            id_select = "SELECT id FROM nonbonded_param WHERE sigma = {0} AND epsilon = {1} " .format(sig, epsi)
            
            cur.execute(nonbonded_param_cmd)
            cur.execute(id_select)
            nbtype = cur.fetchone()[0]
            #print nbtype
            update_particle_charge_nbtype = "UPDATE particle SET nbtype = {0} , charge = {1} WHERE id = {2} " .format (nbtype, chrg, (index-1))
            cur.execute(update_particle_charge_nbtype)
    return [num, charge, sigma, epsilon]
        
#adds new torsion parameters to the dihedral_trig_param table and assigns parameters to the atom combinations in the dihedral_trig_term table
#Also returns the torsion atom1 and atom4 for updating the nonbonded 1-4 interaction tables.
def add_torsion_to_dms(con,paramfile): 
    
    toratm1, toratm2, toratm3, toratm4, tor1, tor2, tor3, tor4 = extract_torsion(paramfile)
    with con:
        cur = con.cursor()
        del_dihedral_param ="DELETE FROM dihedral_trig_param"
        cur.execute(del_dihedral_param)
        for j in range(0,len(toratm1)):
            tatom1 = int(toratm1[j])
            tatom2 = int(toratm2[j])
            tatom3 = int(toratm3[j])
            tatom4 = int(toratm4[j])
            #print tor1[j]
            if float(tor1[j]) == 0.000:
                t1 = float(tor1[j])
            else:
                t1 = float(tor1[j])/2.0
            if float(tor2[j]) == 0.000 :
                t2 = float(tor2[j])
            else:
                t2 = -float(tor2[j])/2.0
            if float(tor3[j]) == 0.000:
                t3 = float(tor3[j])
            else:
                t3 = float(tor3[j])/2.0
            if float(tor4[j]) == 0.000:
                t4 = float(tor4[j])
            else:
                t4 = -float(tor4[j])/2.0

            t0 = float(t1 - t2 + t3 - t4)
            #print t0, t1,t2,t3,t4
            #with con:
            #    cur = con.cursor()
            
            dihedral_param_cmd = """INSERT INTO dihedral_trig_param(fc0, fc1, fc2, fc3, fc4)
            SELECT {0}, {1}, {2}, {3}, {4}
            EXCEPT
            SELECT fc0, fc1, fc2, fc3, fc4 FROM dihedral_trig_param WHERE fc0 = {0} AND fc1 = {1} AND fc2 = {2} AND fc3 = {3} AND fc4 = {4} """ .format(t0, t1, t2, t3, t4)
            
        
            paramid_select = "SELECT id FROM dihedral_trig_param WHERE fc0 = {0} AND fc1 = {1} AND fc2 = {2} AND fc3 = {3} AND fc4 = {4}" .format(t0, t1, t2, t3, t4)
            
            cur.execute(dihedral_param_cmd)
            cur.execute(paramid_select)
            param = cur.fetchone()[0]
            #print param,(tatom1-1),(tatom2-1),(tatom3-1),(tatom4-1)
            update_dihedral_param = "UPDATE dihedral_trig_param SET phi0 = {0} , fc5 = {0}, fc6 = {0} " .format(0.0)
            cur.execute(update_dihedral_param)
            
            update_parent_term = """INSERT INTO dihedral_trig_term(p0, p1, p2, p3)
            SELECT {0}, {1}, {2}, {3}
            EXCEPT
            SELECT p0, p1, p2, p3 FROM dihedral_trig_term WHERE p0 = {0} AND p1 = {1} AND p2 = {2} AND p3 = {3} """ .format((tatom1-1), (tatom2-1), (tatom3-1), (tatom4-1))
            cur.execute(update_parent_term)
            
            update_dihedral_term_param = "UPDATE dihedral_trig_term SET param = {0} WHERE p0 = {1} AND p1 = {2} AND p2 = {3} AND p3 = {4} " .format( param, (tatom1-1), (tatom2-1), (tatom3-1), (tatom4-1))
            cur.execute(update_dihedral_term_param)

            if tatom1 > tatom4 :
                tatom_tmp = tatom4
                tatom4 = tatom1
                tatom1 = tatom_tmp
                
            insert_torsionatm_to_exclusion = """INSERT INTO exclusion(p0,p1)
            SELECT {0}, {1} 
            EXCEPT
            SELECT p0,p1 FROM exclusion WHERE p0 = {0} AND p1 = {1}""" .format(tatom1-1,tatom4-1)

            cur.execute(insert_torsionatm_to_exclusion)
            
    return [toratm1, toratm4]


#updates the nonbonded 1-4 interaction parameters in the pair_12_6_es_param table(calculating aij,bij,qij) and matching id with the atom pairs in the pair_12_6_es_param table.
def add_nb14_to_dms(con, atomtable, toratm1,toratm4, nb14scales_coul, nb14scales_lj):
    num = atomtable[0]
    charge = atomtable[1]
    sigma = atomtable[2]
    epsilon = atomtable[3]
    #toratm1 = torsion[0]
    #toratm4 = torsion[1]
    
    
    #print toratm1, toratm4
    #print len(toratm1)
    #print len(nb14atm1)
    with con:
        cur = con.cursor()
        del_14_es_param = "DELETE FROM pair_12_6_es_param"
        del_14_es_term = "DELETE FROM pair_12_6_es_term"
        cur.execute(del_14_es_term)
        cur.execute(del_14_es_param)
        #number_of_torsions =len(toratm1)
        for i in range(0,len(toratm1)):
            atom1 = int(toratm1[i])
            atom4 = int(toratm4[i])
            chg1 = float(charge[atom1-1])
            sig1 = float(sigma[atom1-1])
            eps1 = float(epsilon[atom1-1])
            chg4 = float(charge[atom4-1])
            sig4 = float(sigma[atom4-1])
            eps4 = float(epsilon[atom4-1])
            q_14 = chg1 * chg4
            sig_14 = math.sqrt(sig1 * sig4)
            eps_14 = math.sqrt(eps1 * eps4)
            #print atom1, atom4
            
            #there are many detected 1,4 interactions in the torsion table which also takes part in 1,3 interactions. So,
            #in order to remove those atom pairs from the 1,4 table, we used the conditional statement where it will only add 
            #the interaction which have a defined 1,4 scale in the nb14 scale table. The others are not added to the table because it
            #already got included in the angle bend table. So, Impact will not calculate them as part of 1,4 interaction.
            try:
                
                s_14_lj = nb14scales_lj[atom1][atom4]
                s_14_coul = nb14scales_coul[atom1][atom4]
                
                #print atom1, atom4, eps1, eps4

                q_14_nb = s_14_coul * q_14
                a_14 = 4 * eps_14 * s_14_lj * (sig_14)**12
                b_14 = 4 * eps_14 * s_14_lj * (sig_14)**6
                
                #print atom1, atom4, s_14_coul, s_14_lj
                es_14_param_cmd = """INSERT INTO pair_12_6_es_param(aij, bij, qij)
                SELECT {0}, {1}, {2}
                EXCEPT
                SELECT aij, bij, qij FROM pair_12_6_es_param WHERE aij = {0} AND bij = {1} AND qij = {2} """ .format(a_14, b_14, q_14_nb)
                
                es_14_idselect = "SELECT id FROM pair_12_6_es_param WHERE aij = {0} AND bij = {1} AND qij = {2} """ .format(a_14, b_14, q_14_nb)
                cur.execute(es_14_param_cmd)
                cur.execute(es_14_idselect)
                idselect_14 = cur.fetchone()[0]
                #print atom1-1, atom4-1, a_14, b_14, idselect_14

                update_es_14_term = """INSERT INTO pair_12_6_es_term(p0, p1, param)
                SELECT {0}, {1}, {2}
                EXCEPT
                SELECT p0, p1, param FROM pair_12_6_es_term WHERE p0 = {0} AND p1 = {1} AND param = {2} """ . format((atom1-1), (atom4-1), idselect_14)
                cur.execute(update_es_14_term)

                #update_es_term_param = "UPDATE pair_12_6_es_term SET param = {0} WHERE p0 = {1} AND p1 = {2} " .format(idselect_14, (atom1-1), (atom4-1))
                cur.execute(update_es_14_term)
                #cur.execute(update_es_term_param)

            except Exception:
                pass
                #print atom1-1, atom4-1
                
            
#for learning purpose: please ignore
"""    print torsion
    print atomtable
    print type(torsion)
    #print len(num)
    nb_14_param = []
    for atom1 in toratm1:
        
        #print tatm1
        for k in range(0, len(num)):
            atomnum = num[k]
            chg = charge[k]
            sig = sigma[k]
            eps = epsilon[k]
            if atom1 == atomnum:
                t1_chg = chg
                #print atom1, t1_chg
                nb_14_param.append([atom1, t1_chg])
    #print nb_14_param

    for atom4 in toratm4:
        for l in range(0, len(num)):
            atomnum = num[l]
            chg = charge[l]
            sig = sigma[l]
            eps = epsilon[l]
            if atom4 == atomnum:
                t4_chg = chg
                nb_14_param.append([atom4, t4_chg])
    print nb_14_param
                
   

for atom1 in toratm1:
        for atom4 in toratm4:
            for k in range(0, len(num)):
                atomnum = num[k]
                chg = charge[k]
                sig = sigma[k]
                eps = epsilon[k]
                if atom1 == atomnum:
                    t1_chg = chg
                if atom4 == atomnum:
                    t4_chg = chg
                    print t4_chg
            #nb_14_param.append((atom1, atom4, t1_chg, t4_chg))
    #print nb_14_param"""
 
            
     
#updates the parameter setting for improper torsions in the dihedral_trig_param and dihedral_trig_term table
def add_improper_torsion_to_dms(con,paramfile):
    im_atom1, im_atom2, im_atom3, im_atom4, im_torsion2 = extract_improper_torsion(paramfile)
    with con:
        cur = con.cursor()
        for m in range(0, len(im_atom1)):
            im_atm1 = int(im_atom1[m])
            im_atm2 = int(im_atom2[m])
            im_atm3 = int(im_atom3[m])
            im_atm4 = int(im_atom4[m])
            if float(im_torsion2[m]) == 0.000:
                im_tor2 = float(im_torsion2[m])
            else:
                im_tor2 = -float(im_torsion2[m])/2.0
            #print im_atm4
            insert_improper_param = """INSERT INTO dihedral_trig_param(fc2)
            SELECT {0}
            EXCEPT 
            SELECT fc2 FROM dihedral_trig_param WHERE fc2 = {0} """ .format(im_tor2)
            
            improper_param_select = "SELECT id FROM dihedral_trig_param WHERE fc2 = {0} " .format(im_tor2)
            cur.execute(insert_improper_param)
            cur.execute(improper_param_select)
            improper_param = cur.fetchone()[0]
            #print improper_param, im_atm1, im_atm2, im_atm3, im_atm4
            update_improper_dihedral_param = "UPDATE dihedral_trig_param SET phi0 = {0}, fc0 = {1}, fc1 = {0}, fc3 = {0}, fc4 = {0} , fc5 = {0} , fc6 = {0} WHERE id = {2} " .format(0.0, -im_tor2, improper_param)
            cur.execute(update_improper_dihedral_param)
            
            update_improper_parent_term = """INSERT INTO dihedral_trig_term(p0, p1, p2, p3)
            SELECT {0}, {1}, {2}, {3}
            EXCEPT
            SELECT p0, p1, p2, p3 FROM dihedral_trig_term WHERE p0 = {0} AND p1 = {1} AND p2 = {2} AND p3 = {3} """ .format((im_atm1-1), (im_atm2-1), (im_atm3-1), (im_atm4-1))
            cur.execute(update_improper_parent_term)

            update_dihedral_term_improper_param = "UPDATE dihedral_trig_term SET param = {0} WHERE p0 = {1} AND p1 = {2} AND p2 = {3} AND p3 = {4} " .format( improper_param, (im_atm1-1), (im_atm2-1), (im_atm3-1), (im_atm4-1))
            cur.execute(update_dihedral_term_improper_param)



#updates the bond stretch parameters into stretch_harm_term and stretch_harm_param tables.
#does not update the 'constrained column in stretch_harm_term table
def add_stretch_to_dms(con,paramfile):
    stretch1, stretch2, k_st, r_eq = extract_stretch(paramfile)
    with con:
        cur = con.cursor()
        del_stretch_param = "DELETE FROM stretch_harm_param"
        cur.execute(del_stretch_param)
        for k in range(0,len(stretch1)):
            st1 = int(stretch1[k])
            st2 = int(stretch2[k])
            k_stretch = float(k_st[k])
            r_eq_stretch = float(r_eq[k])
            
            stretch_param_cmd = """INSERT INTO stretch_harm_param(r0,fc)
            SELECT {0}, {1}
            EXCEPT
            SELECT r0, fc FROM stretch_harm_param WHERE r0 = {0} AND fc = {1} """ .format(r_eq_stretch, k_stretch)

            stretch_param_select = "SELECT id FROM stretch_harm_param WHERE r0 = {0} AND fc = {1} """ .format(r_eq_stretch, k_stretch)
            cur.execute(stretch_param_cmd)
            cur.execute(stretch_param_select)
            stretch_param = cur.fetchone()[0]
            #print stretch_param, st1, st2, k_stretch
            update_parent_stretch_term = """INSERT INTO stretch_harm_term(p0, p1)
            SELECT {0}, {1}
            EXCEPT
            SELECT p0, p1 FROM stretch_harm_term WHERE p0 = {0} AND p1 = {1} """ .format((st1-1), (st2-1))
            cur.execute(update_parent_stretch_term)
            update_stretch_term_param = "UPDATE stretch_harm_term SET param = {0} WHERE p0 = {1} AND p1 = {2} " .format(stretch_param, (st1-1), (st2-1))
            cur.execute(update_stretch_term_param)

            if st1 > st2 :
                st_tmp = st2
                st2 = st1
                st1 = st_tmp
          
            insert_stretch_to_exclusion = """INSERT INTO exclusion(p0,p1)
            SELECT {0}, {1} 
            EXCEPT
            SELECT p0,p1 FROM exclusion WHERE p0 = {0} AND p1 = {1}""" .format(st1-1,st2-1)

            cur.execute(insert_stretch_to_exclusion)
            


#updates the bend angle parameters in angle_harm_term and angle_harm_param tables.
def add_bend_to_dms(con,paramfile):
    bendatm1, bendatm2, bendatm3, bend_k, bend_theta_eq = extract_bend(paramfile)
    with con:
        cur = con.cursor()
        del_angle_param = "DELETE FROM angle_harm_param"
        cur.execute(del_angle_param)
        for l in range(0,len(bendatm1)):
            batm1 = int(bendatm1[l])
            batm2 = int(bendatm2[l])
            batm3 = int(bendatm3[l])
            k_bend = float(bend_k[l])
            theta_eq_bend = float(bend_theta_eq[l])
            
            angle_param_cmd = """INSERT INTO angle_harm_param(theta0, fc)
            SELECT {0}, {1}
            EXCEPT
            SELECT theta0, fc FROM angle_harm_param WHERE theta0 = {0} AND fc = {1} """ .format(theta_eq_bend, k_bend)
            angle_param_select = "SELECT id FROM angle_harm_param WHERE theta0 = {0} AND fc = {1} " .format(theta_eq_bend, k_bend)
            cur.execute(angle_param_cmd)
            cur.execute(angle_param_select)
            angle_param = cur.fetchone()[0]
            #print angle_param, batm1,batm2,batm3
            
            update_parent_angle_term = """INSERT INTO angle_harm_term(p0, p1, p2)
            SELECT {0}, {1}, {2}
            EXCEPT
            SELECT p0, p1, p2 FROM angle_harm_term WHERE p0 ={0} AND p1 = {1} AND p2 = {2} """. format((batm1-1), (batm2-1), (batm3-1))
            cur.execute(update_parent_angle_term)
            
            update_angle_term_param = "UPDATE angle_harm_term SET param = {0} WHERE p0 = {1} AND p1 = {2} AND p2 = {3} """ .format(angle_param, (batm1-1), (batm2-1), (batm3-1))
            cur.execute(update_angle_term_param)

            if batm1 > batm3 :
                batm_tmp = batm3
                batm3 = batm1
                batm1 = batm_tmp

            insert_angle_to_exclusion = """INSERT INTO exclusion(p0,p1)
            SELECT {0}, {1} 
            EXCEPT
            SELECT p0,p1 FROM exclusion WHERE p0 = {0} AND p1 = {1}""" .format(batm1-1,batm3-1)

            cur.execute(insert_angle_to_exclusion)
    

#main function which reads the parameter and calls different functions to parse and update the database file
def update_dms(param, jobname):
    opls2005_dms = jobname + ".dms"
    dms = jobname + "_opls3.dms"
    #dms file is copied at the beginning to keep the original opls2005 parameters for comparison
    shutil.copyfile(opls2005_dms, dms)
    #sql connection opened with the copied dms file
    con = lite.connect(dms)
    with con:
        cur = con.cursor()
        del_exclusion = "DELETE FROM exclusion"
        cur.execute(del_exclusion)
        
    #read parameter file    
    with open(param, 'r') as infile:
        
        for line in infile:
            if line.startswith('Title'):
                pass
            row = line.split()
            if not row:
                pass
                
            elif row[0] == 'Atom':
                num, charge, sigma, epsilon = extract_atomtable(infile)
                atom_table = add_atomtable_to_dms(con, num, charge, sigma, epsilon)
            
            elif row[0] == 'Stretch':
                add_stretch_to_dms(con,infile)
                
            elif row[0] == 'Bend_angle':
                add_bend_to_dms(con,infile)
                
            elif row[0] == 'Proper_torsion':
                torsionatm1, torsionatm4 = add_torsion_to_dms(con,infile)
                
            
            elif row[0] == 'Improper_torsion':
                add_improper_torsion_to_dms(con,infile)

            elif row[0] == 'NB14_pair':
                nb14scales_coul, nb14scales_lj = extract_nb14_scaling_factor(num,infile)
                add_nb14_to_dms(con, atom_table, torsionatm1, torsionatm4, nb14scales_coul, nb14scales_lj)
                #del_14(con)
            
                
def main():
    try:
        param = sys.argv[1]
        jobname = sys.argv[2]
    except IndexError as e:
        print "Usage: python <file.py> <paramfile> <jobname without file extension>"
        print "Please provide all the arguments... exiting.."
        print e
        sys.exit(1)
    
    update_dms(param, jobname)


if __name__ == "__main__":
    main()
    
