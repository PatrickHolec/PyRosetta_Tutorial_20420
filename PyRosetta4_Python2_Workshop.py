# edits by Corban Swain - 10/3/2017
# updating 20.420 PyRosetta tutorial for Pyrosetta 4 and python 2.7
# original tutorial is here
# https://github.com/PatrickHolec/Pyrosetta_Tutorial_20420

#### to run this script ####
# 1. Download and install Python 2.7
# 2. Download and extract PyRosetta4 for python 2.7
# 3. in Terminal, cd to wherever PyRosetta is extracted
# 4. in Terminal run the following commands to install PyRosetta4
#    --> $ cd setup
#    --> $ sudo python2 setup.py install
# 5. Navigate to https://github.com/PatrickHolec/Pyrosetta_Tutorial_20420
#    and download the repository .... OR clone the repo using git via terminal
# 6. Copy this script file into the the same folder as the tutorial repo and run it
#    using python 2

# This code is UNCOMMENTED ... except for the few places I noticed
# something.  To see explainations of the steps, reference the TA's
# Jypyter tutorial.  All the function calls are in the same order.

# most of the changes are editing where the functions and classes are
# located there are a couple of syntatical changes.  I have NOT
# highlighted any of these changes here...you'll have to visually look
# at this and the previous Jupyter tutorial to see the difference.

# Also, I'm not one of the cool kids and unfortunately didnt know how
# to use Jupyter, (go emacs!) so this is written as a bunch of
# functions in a .py file.

# -CS


from pyrosetta import *

def initialize():
    init()

def load():
    import subprocess
    subprocess.call('babel pyrrolysine.pdb pyrrolysine.mdl', shell=True)
    subprocess.call('molfile2params/molfile_to_params.py pyrrolysine.mdl', shell=True)

def importPdb():
    global pose, scorefxn
    params_list = ['LG.params']
    pose = Pose()
    generate_nonstandard_residue_set(pose, params_list)
    pose_from_file(pose,'2AVX_pyrrolysine.pdb')
    print pose
    scorefxn = get_fa_scorefxn()



def printPyRoseDir():
    import os
    print os.path.dirname(rosetta.__file__)
    # this prints: /Library/Python/2.7/site-packages/pyrosetta-4.0-py2.7-macosx-10.13-intel.egg/pyrosetta
    # ^ dont use this
    # ^ this is where python has loaded in PyRosetta
    #
    # instead use: /Users/corbann.swain/Applications/PyRosetta4.Release.python27.mac.release-1
    # ^ this is wherever you put the PyRosetta4 folder
    #
    # In PyMOL:
    #
    # PyMOL> cd /Users/corbann.swain/Applications/PyRosetta4.Release.python27.mac.release-154
    # PyMOL> run PyMOL-RosettaServer.py
    #
    # Done in PyMOL.

def sendToPyMOL():
    global pymol
    pymol = PyMOLMover()
    pymol.apply(pose)

    print 'My score: %s' % scorefxn(pose)
    pymol.send_energy(pose)

def showHBonds():
    pymol.send_hbonds(pose)

def findMimima():
    min_mover = rosetta.protocols.simple_moves.MinMover()
    move_map = MoveMap()
    move_map.set_bb_true_range(40,130)
    min_mover.movemap(move_map)
    min_mover.score_function(scorefxn)

    pose_min_move = Pose()
    pose_min_move.assign(pose)

    observer = rosetta.protocols.moves.AddPyMOLObserver(pose, True)

    min_mover.apply(pose_min_move)

    print 'Original Energy: %.2f' % scorefxn(pose)
    print 'Energy after min_mover: %.2f' % scorefxn(pose_min_move)
    
def fastRelax():
    global fast_relax
    fast_relax =rosetta.protocols.relax.FastRelax()
    fast_relax.set_scorefxn(scorefxn)

    pose_fast_relax = Pose()
    pose_fast_relax.assign(pose)

    pose_fast_relax.pdb_info().name('fast_relax')
    print pose.pdb_info().name
    print pose_fast_relax.pdb_info().name()

    observer =rosetta.protocols.moves.AddPyMOLObserver(pose_fast_relax, True)

    fast_relax.apply(pose_fast_relax)

    print 'Finished!'

def fastRelax2():
    pose_fast_relax = Pose()
    pose_fast_relax.assign(pose)
    fast_relax.apply(pose_fast_relax)

    print 'Original Energy: %.2f' % scorefxn(pose)
    print 'Energy after fast relax: %.2f' % scorefxn(pose_fast_relax)

    pose_fast_relax.pdb_info().name('fast_relax_final')
    pymol.apply(pose_fast_relax)
    pymol.send_energy(pose_fast_relax)
    
def mutataRes():
    pose_pack_mover = Pose()
    pose_pack_mover.assign(pose)
    pose_pack_mover.pdb_info().name('pack_mover')

    task = standard_packer_task(pose_pack_mover)
    task.temporarily_fix_everything()
    task.temporarily_set_pack_residue(124, True)
    task.temporarily_set_pack_residue(123, True)

    pack_mover = rosetta.protocols.simple_moves.PackRotamersMover(scorefxn, task)

    observer = rosetta.protocols.moves.AddPyMOLObserver(pose_pack_mover, True)
    pack_mover.apply(pose_pack_mover)

    print 'Original Score: %s' % scorefxn(pose)
    print 'Score after mutating 124, 123: %s' % scorefxn(pose_pack_mover)

def controlMut():
    from pyrosetta.toolbox import generate_resfile_from_pdb, generate_resfile_from_pose
    # generate_resfile_from_pdb('2AVX_pyrrolysine.pdb', 'my.resfile')
    #  ^ causes error  EXCN_Base::what()
    # likely due to uncleaned file?
    #
    # what if we clean first ...
    from pyrosetta.toolbox import cleanATOM
    cleanATOM('2AVX_pyrrolysine.pdb')
    generate_resfile_from_pdb('2AVX_pyrrolysine.clean.pdb','my2.resfile')
    # ^ doesnt work either causes error: ERROR: Unrecognized residue: LG1

    # The following does work.
    generate_resfile_from_pose(pose, 'my.resfile')

def editResfile():
    import re
    global new_rfile
    fname = 'my'
    fext = 'resfile'
    rfile = open('.'.join([fname, fext]))
    text = rfile.read()
    rfile.close()
    res_to_mutate = [89, 101, 123, 139, 147]
    reg = r'( *(%s)\s*A\s*)(\w+)( *)' % ' | '.join(map(str,res_to_mutate))
    subStr = r'\1PIKAA STNQ\4'
    def subFun(m): return m.expand(subStr)         
    new_text = re.sub(reg,subFun,text)
    new_rfile = open('.'.join([fname + '_MUT', fext]),'w+')
    new_rfile.write(new_text)
    new_rfile.close()

def packMutant():
    pose_pack_mover = Pose()
    pose_pack_mover.assign(pose)
    pose_pack_mover.pdb_info().name('defined_pack_mover')

    observer = rosetta.protocols.moves.AddPyMOLObserver(pose_pack_mover, True)

    kT = 1.0
    mc = MonteCarlo(pose_pack_mover, scorefxn, kT)

    for i in range(11):
        task = rosetta.core.pack.task.TaskFactory.create_packer_task(pose_pack_mover)
        rosetta.core.pack.task.parse_resfile(pose_pack_mover, task, new_rfile.name)
        packer_task = rosetta.protocols.simple_moves.PackRotamersMover(scorefxn, task)

        min_mover = rosetta.protocols.simple_moves.MinMover()
        mm70150 = MoveMap()
        mm70150.set_bb_true_range(70, 150)
        min_mover.movemap(mm70150)
        min_mover.score_function(scorefxn)

        seq_mover = SequenceMover()
        seq_mover.add_mover(packer_task)
        seq_mover.add_mover(min_mover)

        trial_pack_min_mover = TrialMover(seq_mover, mc)

        trial_pack_min_mover.apply(pose_pack_mover)
        print 'Score: %s' % scorefxn(pose_pack_mover)

    # PyMOL> select my_ser, (resi 89+101+123+139+147) & defined_pack_mover

def repeatMovers():
    pose_repeat_pack_mover = Pose()
    pose_repeat_pack_mover.assign(pose)
    pose_repeat_pack_mover.pdb_info().name('repeat_packer')

    kT = 1.0
    movemap = MoveMap()
    movemap.set_bb(True)
    small_mover = rosetta.protocols.simple_moves.SmallMover(movemap, kT, 1)

    mc = MonteCarlo(pose_repeat_pack_mover, scorefxn, kT)
    trial_mover = TrialMover(small_mover, mc)

    n = 20
    repeat_mover = RepeatMover(trial_mover, n)

    repeat_mover.apply(pose_repeat_pack_mover)

    print 'original: %s' % scorefxn(pose)
    print 'new: %s' % scorefxn(pose_repeat_pack_mover)
    # ^ all of this code computes, but it doesnt seem to do anything...not sure whats wrong.

def main():
    # comment out sections you dont want to run...  for a given python
    # session, any function will run after all the previous functions
    # have run once
    
    initialize()
    load()
    importPdb()
    printPyRoseDir()
    sendToPyMOL()
    showHBonds()
    findMimima()
    
    # Both of the fast relax sessions take a while...uncomment to run
    # fastRelax()
    # fastRelax2()
    
    mutataRes()
    controlMut()
    editResfile()
    packMutant()
    repeatMovers()

main()
    
