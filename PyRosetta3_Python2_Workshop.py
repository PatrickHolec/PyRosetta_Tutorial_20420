
"""

Athena cluster - PyRosetta Script

"""


print '-'*40

print 'Welcome to the PyRosetta tutorial for 20.420!'
print 'Since we cannot get users to remotely connect to Jupyter, we hope this program will help demonstrate PyRosetta usage without this dynamic environment.'
print 'First- PyRosetta is going to load the PyRosetta library which may take sometime (~1 minute)'
print 'Do not freak out if it seems like this script is frozen!'

from rosetta import *
init()

print 'PyRosetta library loaded!'
print '-'*40

cmd_descriptions = {
        'cmd1':'load PDB w/o adding ligand params',
        'cmd2':'convert ligand into a mdl file',
        'cmd3':'load PDB w/ adding ligand params',
        'cmd4':'N/A',
        'cmd5':'N/A',
        'cmd6':'N/A',
        'cmd7':'Find energy of structure after repacking structure (simple)', 'cmd8':'Create a fast relax object (w/o PyMOL link)',
        'cmd9':'Print out PyRosetta scores on fast relax/non-relax structures',
        'cmd10':'N/A',
        'cmd11':'Mutate residues 123/124 and repack structure',
        'cmd12':'Generate resfiles for structure',
        'cmd13':'Use my.resfile to mutate protein and find local minimum in energy',
        'cmd14':'Use repeat movers to simulate mutations on a loop'
            }

cmd_all = ['cmd'+str(i) for i in xrange(1,len(cmd_descriptions)+1)]


print 'To use this script, just enter commands indicated in the Jupyter notebook. These are:'
for k in cmd_all:
    print '> {} - {}'.format(k,cmd_descriptions[k])
print '-'*40


"""
Start actually command line prompt
"""

while True:

    cmd = raw_input('> ')
   
    try: # catch for KeyboardInterupt!

        ###########################
        if cmd == 'cmd1':

            try:
                print 'Running Python commands:'
                print "> pose = pose_from_pdb('2AVX_pyrrolysine.pdb')"
                print '-'*40

                pose = pose_from_pdb('2AVX_pyrrolysine.pdb')

                print '-'*15 + ' FINISHED COMMAND ' + '-'*15
            except RuntimeError:
                print 'We ran into an error (on purpose, though)!'
                print '-'*15 + ' FINISHED COMMAND ' + '-'*15

        ###########################
        elif cmd == 'cmd2':
            
            print 'Running Python commands:'
            print "> import subprocess"
            print "> subprocess.call('babel pyrrolysine.pdb pyrrolysine.mdl', shell=True)"
            print "> subprocess.call('molfile2params/molfile_to_params.py pyrrolysine.mdl',shell=True)"
            print '-'*40

            import subprocess
            subprocess.call('babel pyrrolysine.pdb pyrrolysine.mdl', shell=True)
            subprocess.call('molfile2params/molfile_to_params.py pyrrolysine.mdl',shell=True)

            print '-'*15 + ' FINISHED COMMAND ' + '-'*15

        ###########################
        elif cmd == 'cmd3':
            
            print 'Running Python commands:'
            print "> # Load up the parameters"
            print "> params_list = Vector1(['LG.params'])"
            print "> res_set = generate_nonstandard_residue_set(params_list)"
            print "> # Import your ligand-receptor pdb file"
            print "> pose = pose_from_pdb('2AVX_pyrrolysine.pdb')"
            print "> scorefxn = get_fa_scorefxn() # The standard full-atom scorefunction"

            print '-'*40

            params_list = Vector1(['LG.params'])
            res_set = generate_nonstandard_residue_set(params_list)
            pose = pose_from_pdb('2AVX_pyrrolysine.pdb')
            scorefxn = get_fa_scorefxn() # The standard full-atom scorefunction

            print '-'*15 + ' FINISHED COMMAND ' + '-'*15
        
        ###########################
        elif cmd == 'cmd7':
            
            print 'Running Python commands:'

            print "> # Let's find a local minima"
            print "> min_mover = MinMover()"
            print "> # Set which residues are flexible with a move map"
            print "> move_map = MoveMap()"
            print "> move_map.set_bb_true_range(40,130)"
            print "> # Give your movemap and scorefunction to the Mover"
            print "> min_mover.movemap(move_map)"
            print "> min_mover.score_function(scorefxn)"

            print "> # Keep track of what's happening to visualize in pymol"
            print "> pose_min_move = Pose()"
            print "> pose_min_move.assign(pose)"

            print "> # This should keepy send a frame to pymol every time the pose changes"
            print "> # so you can actually see what's going on. Caution that I don't know"
            print "> # if it is only keeping the poses accepted by the Metropolis criteria,"
            print "> # or if it keeps any pose at all. (I think it's keeping anything"
            print "> # at all)"

            print "> min_mover.apply(pose_min_move)"

            print "> print 'Original Energy: %.2f' % scorefxn(pose)"
            print "> print 'Energy after min_mover: %.2f' % scorefxn(pose_min_move)"

            print '-'*40

            # Let's find a local minima
            min_mover = MinMover()
            # Set which residues are flexible with a move map
            move_map = MoveMap()
            move_map.set_bb_true_range(40,130)
            # Give your movemap and scorefunction to the Mover
            min_mover.movemap(move_map)
            min_mover.score_function(scorefxn)

            # Keep track of what's happening to visualize in pymol
            pose_min_move = Pose()
            pose_min_move.assign(pose)

            min_mover.apply(pose_min_move)

            print 'Original Energy: %.2f' % scorefxn(pose)
            print 'Energy after min_mover: %.2f' % scorefxn(pose_min_move)

            print '-'*15 + ' FINISHED COMMAND ' + '-'*15


        ###########################
        elif cmd == 'cmd8':
            
            print 'Running Python commands:'

            print "> fast_relax = FastRelax()"
            print "> fast_relax.set_scorefxn(scorefxn)"
            print "> # Again, make a pose and follow it with pymol_observer"
            print "> pose_fast_relax = Pose()"
            print "> pose_fast_relax.assign(pose)"
            print "> # rename it so that PyMol won't overwrite our "
            print "> pose_fast_relax.pdb_info().name('fast_relax')"
            print "> print pose.pdb_info().name()"
            print "> print pose_fast_relax.pdb_info().name()"

            print "> # Relax"
            print "> fast_relax.apply(pose_fast_relax)"

            print "> print 'Finished!'"

            print '-'*40

            fast_relax = FastRelax()
            fast_relax.set_scorefxn(scorefxn)
            # Again, make a pose and follow it with pymol_observer
            pose_fast_relax = Pose()
            pose_fast_relax.assign(pose)
            # rename it so that PyMol won't overwrite our 
            pose_fast_relax.pdb_info().name('fast_relax')
            print pose.pdb_info().name()
            print pose_fast_relax.pdb_info().name()

            # Relax
            fast_relax.apply(pose_fast_relax)

            print 'Finished!'


            print '-'*15 + ' FINISHED COMMAND ' + '-'*15


        ###########################
        elif cmd == 'cmd9':
            
            print 'Running Python commands:'

            print "> # Run FastRelax by itself"
            print "> pose_fast_relax = Pose()"
            print "> pose_fast_relax.assign(pose)"
            print "> fast_relax.apply(pose_fast_relax)"

            print "> print 'Original Energy: %.2f' % scorefxn(pose)"
            print "> print 'Energy after fast relax: %.2f' % scorefxn(pose_fast_relax)"

            print '-'*40

            # Run FastRelax by itself
            pose_fast_relax = Pose()
            pose_fast_relax.assign(pose)
            fast_relax.apply(pose_fast_relax)

            print 'Original Energy: %.2f' % scorefxn(pose)
            print 'Energy after fast relax: %.2f' % scorefxn(pose_fast_relax)

            print '-'*15 + ' FINISHED COMMAND ' + '-'*15


        ###########################
        elif cmd == 'cmd11':
            
            print 'Running Python commands:'

            print "> # Copy starting pose"
            print "> pose_pack_mover = Pose()"
            print "> pose_pack_mover.assign(pose)"
            print "> pose_pack_mover.pdb_info().name('pack_mover')"

            print "> # Make a Pack Mover"
            print "> # Default is that all residues can mutate and repack. "
            print "> # Restrict to only Residues 124, 123"
            print "> task = standard_packer_task(pose_pack_mover)"
            print "> task.temporarily_fix_everything()"
            print "> task.temporarily_set_pack_residue(124, True)"
            print "> task.temporarily_set_pack_residue(123, True)"

            print "> # Make the mover with your restrictions "
            print "> pack_mover = PackRotamersMover(scorefxn, task)"
            print "> # Repack"
            print "> pack_mover.apply(pose_pack_mover)"

            print "> print 'Original Score: %s' % scorefxn(pose)"
            print "> print 'Score after mutating 124, 123: %s' % scorefxn(pose_pack_mover)"

            print '-'*40

            # Copy starting pose
            pose_pack_mover = Pose()
            pose_pack_mover.assign(pose)
            pose_pack_mover.pdb_info().name('pack_mover')

            # Make a Pack Mover
            # Default is that all residues can mutate and repack. 
            # Restrict to only Residues 124, 123
            task = standard_packer_task(pose_pack_mover)
            task.temporarily_fix_everything()
            task.temporarily_set_pack_residue(124, True)
            task.temporarily_set_pack_residue(123, True)

            # Make the mover with your restrictions 
            pack_mover = PackRotamersMover(scorefxn, task)
            # Repack
            pack_mover.apply(pose_pack_mover)

            print 'Original Score: %s' % scorefxn(pose)
            print 'Score after mutating 124, 123: %s' % scorefxn(pose_pack_mover)

            print '-'*15 + ' FINISHED COMMAND ' + '-'*15

        ###########################
        elif cmd == 'cmd12':
            
            print 'Running Python commands:'

            print "> # Generate a resfile from a pdb or a pose"
            print "> from toolbox import generate_resfile_from_pdb, generate_resfile_from_pose"
            print "> generate_resfile_from_pdb('2AVX_pyrrolysine.pdb', 'my.resfile')"
            print "> generate_resfile_from_pose(pose, 'my.resfile')"

            print '-'*40

            # Generate a resfile from a pdb or a pose
            from toolbox import generate_resfile_from_pdb, generate_resfile_from_pose
            generate_resfile_from_pdb('2AVX_pyrrolysine.pdb', 'my.resfile')
            generate_resfile_from_pose(pose, 'my.resfile')

            print '-'*15 + ' FINISHED COMMAND ' + '-'*15

        ###########################
        elif cmd == 'cmd13':
            
            print 'Running Python commands:'

            print "> # Go and edit the resfile"

            print "> # Copy starting pose"
            print "> pose_pack_mover = Pose()"
            print "> pose_pack_mover.assign(pose)"
            print "> pose_pack_mover.pdb_info().name('defined_pack_mover')"

            print "> # monte carlo object - will allow a TrialMover to "
            print "> # accept or reject moves based on metropolis criteria"
            print "> kT = 1.0"
            print "> mc = MonteCarlo(pose_pack_mover, scorefxn, kT)"

            print "> # Run ten times"
            print "> for i in range(0,11):"
            print ">     # Make the Pack Mover by creating the task (which AA's can be"
            print ">     # mutated), and giving it to the mover"
            print ">     task = TaskFactory.create_packer_task(pose_pack_mover)"
            print ">     parse_resfile(pose_pack_mover, task, 'my.resfile')"
            print ">     packer_task = PackRotamersMover(scorefxn, task)"
                
            print ">     # Make minimization mover"
            print ">     min_mover = MinMover()"
            print ">     mm70150 = MoveMap()"
            print ">     mm70150.set_bb_true_range(70, 150)"
            print ">     min_mover.movemap(mm70150)"
            print ">     min_mover.score_function(scorefxn)"
                
            print ">     # String the packer and minimization movers together"
            print ">     seq_mover = SequenceMover()"
            print ">     seq_mover.add_mover(packer_task)"
            print ">     seq_mover.add_mover(min_mover)"

            print ">     # And make them satisfy a metropolis montecarlo category"
            print ">     trial_pack_min_mover = TrialMover(seq_mover, mc)"
                
            print ">     #profanity"
            print ">     trial_pack_min_mover.apply(pose_pack_mover)"
            print ">     print 'Score: %s' % scorefxn(pose_pack_mover)"

            print '-'*40

            # Go and edit the resfile

            # Copy starting pose
            pose_pack_mover = Pose()
            pose_pack_mover.assign(pose)
            pose_pack_mover.pdb_info().name('defined_pack_mover')

            # monte carlo object - will allow a TrialMover to 
            # accept or reject moves based on metropolis criteria
            kT = 1.0
            mc = MonteCarlo(pose_pack_mover, scorefxn, kT)

            # Run ten times
            for i in range(0,11):
                # Make the Pack Mover by creating the task (which AA's can be
                # mutated), and giving it to the mover
                task = TaskFactory.create_packer_task(pose_pack_mover)
                parse_resfile(pose_pack_mover, task, 'my.resfile')
                packer_task = PackRotamersMover(scorefxn, task)
                
                # Make minimization mover
                min_mover = MinMover()
                mm70150 = MoveMap()
                mm70150.set_bb_true_range(70, 150)
                min_mover.movemap(mm70150)
                min_mover.score_function(scorefxn)
                
                # String the packer and minimization movers together
                seq_mover = SequenceMover()
                seq_mover.add_mover(packer_task)
                seq_mover.add_mover(min_mover)
                
                # And make them satisfy a metropolis montecarlo category
                trial_pack_min_mover = TrialMover(seq_mover, mc)
                
                #profanity
                trial_pack_min_mover.apply(pose_pack_mover)
                print 'Score: %s' % scorefxn(pose_pack_mover)

            print '-'*15 + ' FINISHED COMMAND ' + '-'*15

        ###########################
        elif cmd == 'cmd14':
            
            print 'Running Python commands:'

            print "> # Repeat Mover Example"

            print "> pose_repeat_pack_mover = Pose()"
            print "> pose_repeat_pack_mover.assign(pose)"
            print "> pose_repeat_pack_mover.pdb_info().name('repeat_packer')"

            print "> # Take a backbone mover"
            print "> kT = 1.0"
            print "> movemap = MoveMap()"
            print "> movemap.set_bb(True)"
            print "> small_mover = SmallMover(movemap, kT, 1)"

            print "> # Combine it witha montecarlo object, so that only moves that "
            print "> # satisfy the metropolis criteria are kept"
            print "> mc = MonteCarlo(pose_repeat_pack_mover, scorefxn, kT)"
            print "> trial_mover = TrialMover(small_mover, mc)"

            print "> # Then repeat this move 20 times"
            print "> n=20"
            print "> repeat_mover = RepeatMover(trial_mover, n)"
            print "> repeat_mover.apply(pose_repeat_pack_mover)"
            print "> print 'original: %s' % scorefxn(pose)"
            print "> print 'new: %s' % scorefxn(pose_repeat_pack_mover)"


            print '-'*40

            # Repeat Mover Example

            pose_repeat_pack_mover = Pose()
            pose_repeat_pack_mover.assign(pose)
            pose_repeat_pack_mover.pdb_info().name('repeat_packer')

            # Take a backbone mover
            kT = 1.0
            movemap = MoveMap()
            movemap.set_bb(True)
            small_mover = SmallMover(movemap, kT, 1)

            # Combine it witha montecarlo object, so that only moves that 
            # satisfy the metropolis criteria are kept
            mc = MonteCarlo(pose_repeat_pack_mover, scorefxn, kT)
            trial_mover = TrialMover(small_mover, mc)

            # Then repeat this move 20 times
            n=20
            repeat_mover = RepeatMover(trial_mover, n)
            repeat_mover.apply(pose_repeat_pack_mover)
            print 'original: %s' % scorefxn(pose)
            print 'new: %s' % scorefxn(pose_repeat_pack_mover)


            print '-'*15 + ' FINISHED COMMAND ' + '-'*15

        elif cmd in ['cmd4','cmd5','cmd6','cmd10']:

            print '-'*40
            print 'Warning: this code cannot be executed remotely in Athena!'
            print 'Skipping...'
            print '-'*15 + ' FINISHED COMMAND ' + '-'*15

        elif cmd.lower() == 'exit':
            print 'Exiting PyRosetta tutorial...'
            break

        elif cmd == '':
            continue # no user input

        else:
            print 'Command unrecognized ({})'.format(cmd)
            print 'Returning to prompt...'
            print '-'*15 + ' FINISHED COMMAND ' + '-'*15

    except KeyboardInterrupt:
        print 'User interrupted program! Returning to main...'
        print '-'*40


