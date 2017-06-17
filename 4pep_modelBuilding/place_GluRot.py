import sys, os, numpy as np 
from prody import * 


# Grab rotamer pdbs from directory 
rotDir 		= sys.argv[2]
rot4 		= parsePDB( os.path.join( rotDir, 'glu_rot4h.pdb' ) ).select( 'not hydrogen' )
# rot 13, 6 
rot5 		= parsePDB( os.path.join( rotDir, 'glu_rot5h.pdb' ) ).select( 'not hydrogen' )
rot3 		= parsePDB( os.path.join( rotDir, 'glu_rot3h.pdb' ) ).select( 'not hydrogen' )


rotDict 	= { "A" : [ rot4 ], "B" : [ rot3, rot5 ] } 

# Grab tetra-acetate diRh pdb
diRh 		= parsePDB( sys.argv[3] ) 



# Hash which rotamers to try at each chain... 
residue  = 12
rot_dict = {}

# Look through directory of 
bbDir = sys.argv[1]
for i in os.listdir( bbDir ): 
	if i[-3:] != 'pdb': continue

	bb_Path = os.path.join( bbDir, i )
	pdb 	= parsePDB( bb_Path ).select( 'not hydrogen' )

	### chain A Glu
	mobile 	= rot4.select( 'bb' ).copy()
	newRot  = rot4.copy()
	target  = pdb.select( 'bb chain A not hydrogen resnum 12'  ) 
	newRot.setChids( [ 'A' for x in newRot.getChids()] )
	newRot.setResnums( [ '12' for x in newRot.getResnums()] )
	newRot.setSegnames( [ 'A' for x in newRot.getSegnames()] )
	CD_a 	= newRot.select( 'name CD').copy()

	trans	= calcTransformation( mobile, target ) 
	newGluA = applyTransformation( trans, newRot )

	newChainA = pdb.select( 'chain A resnum 1 to 11' ).copy() + newRot + pdb.select( 'chain A resnum 13 to 25' ).copy() 
	#####


	### chain B Glu
	mobile 	= rot3.select( 'bb' ).copy()
	newRot  = rot3.copy()
	target  = pdb.select( 'bb chain B not hydrogen resnum 12'  ) 

	trans	= calcTransformation( mobile, target ) 
	newGluA = applyTransformation( trans, newRot )
	newRot.setChids( [ 'B' for x in newRot.getChids()] )
	newRot.setResnums( [ '12' for x in newRot.getResnums()] )
	newRot.setSegnames( [ 'B' for x in newRot.getSegnames()] )
	CD_b 	= newRot.select( 'name CD').copy()

	newChainB = pdb.select( 'chain B resnum 1 to 11' ).copy() + newRot + pdb.select( 'chain B resnum 13 to 25' ).copy() 
	#####


	# Quit if CG residues are too far away from each other. 
	dist =  calcDistance( CD_a, CD_b )
	if dist > 7:
		continue
	else: 
		print bb_Path, round( dist[0], 1 )


	### chain C Glu
	mobile 	= rot4.select( 'bb' ).copy()
	newRot  = rot4.copy()
	target  = pdb.select( 'bb chain C not hydrogen resnum 12'  ) 

	trans	= calcTransformation( mobile, target ) 
	newGluA = applyTransformation( trans, newRot )
	newRot.setChids( [ 'C' for x in newRot.getChids()] )
	newRot.setResnums( [ '12' for x in newRot.getResnums()] )
	newRot.setSegnames( [ 'C' for x in newRot.getSegnames()] )


	newChainC = pdb.select( 'chain C resnum 1 to 11' ).copy() + newRot + pdb.select( 'chain C resnum 13 to 25' ).copy() 
	#####


	### chain D Glu
	mobile 	= rot4.select( 'bb' ).copy()
	newRot  = rot4.copy()
	target  = pdb.select( 'bb chain D not hydrogen resnum 12'  ) 

	trans	= calcTransformation( mobile, target ) 
	newGluA = applyTransformation( trans, newRot )
	newRot.setChids( [ 'D' for x in newRot.getChids()] )
	newRot.setResnums( [ '12' for x in newRot.getResnums()] )
	newRot.setSegnames( [ 'D' for x in newRot.getSegnames()] )


	newChainD = pdb.select( 'chain D resnum 1 to 11' ).copy() + newRot + pdb.select( 'chain D resnum 13 to 25' ).copy() 
	#####





	newPdb = newChainA + newChainB + newChainC + newChainD

#	writePDB( 'tmpAB.pdb', newPdb )


#	sys.exit()