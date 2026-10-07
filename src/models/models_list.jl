##########################################################################################################################################
############################################################### MODELS ##################################################################
##########################################################################################################################################

# Curie Weiss with external field 
include("Ferro/CurieWeiss.jl")
include("Ferro/VectorCW.jl")


# coupled CurieWeiss
include("Ferro/RepCW.jl")
include("Ferro/RepCWHomogeneous.jl")
include("Ferro/RepCWHomogeneous_betac.jl")


# RFIM
include("RFIM.jl")

# SK with external field and coupling mean
include("SK.jl")

# Bethe approximation in ferromagnetic models on lattices
include("Ferro/BetheFerro.jl")



################ Hopfield  ###############

include("Hopfield/HopfieldSignal_utils.jl")
# Standard Hopfield low load regime with mixed states retrieval (arbitrary l ) -> independent and binary patterns for now
include("Hopfield/HopfieldLL.jl")
include("Hopfield/HopfieldHL.jl")


# Standard Hopfield in the high load with pure states retrieval only (l=1)

include("Hopfield/HopfieldHLT0.jl")
include("Hopfield/HopfieldHL1RSB.jl")
include("Hopfield/HopfieldHL1RSBT0.jl")

# Hopfield with Gaussian noise -> low load (equivalent to SK when retrieving one pattern)
include("Hopfield/HopfieldLL_GaussNoise.jl")



################ BAM ###############
include("BAM/BAM_utils.jl")
include("BAM/BAMLL.jl")
include("BAM/BAMHL.jl")
include("BAM/BAMHLT0.jl")
include("BAM/BAMHLT0RSstandard.jl")
include("BAM/BAMHLT0Mixtures.jl")
include("BAM/BAMHL1RSBT0.jl")
include("BAM/HopfieldTauBAMK2.jl")

################### Replicated Hopfields #########################



############################# RBMs ############################
include("RestrictedBoltzmannMachines/RBMBinaryBinaryLL_signalutils.jl")
include("RestrictedBoltzmannMachines/RBMBinaryBinaryLL.jl")
include("RestrictedBoltzmannMachines/RBMBinaryBinaryLL_onlyPure.jl")


include("RestrictedBoltzmannMachines/RBMBinaryReLU.jl")
include("RestrictedBoltzmannMachines/RBMBinaryReLUT0.jl")
include("RestrictedBoltzmannMachines/RBMBinaryReLUDiluted.jl")
include("RestrictedBoltzmannMachines/RBMBinaryReLUDilutedexpanded.jl")
include("RestrictedBoltzmannMachines/RBMLowLoadRepeatedHidden.jl")


######################################## PERCEPTRON WITH POTENTIALS  #####################################
abstract type Perceptron <: FPModel end

abstract type Potential end

include("Perceptron_with_Potentials/Potentials/potentials.jl")
include("Perceptron_with_Potentials/perceptron_utils.jl")

include("Perceptron_with_Potentials/BinaryPerceptronRS.jl")
include("Perceptron_with_Potentials/BinaryPerceptron1RSB.jl")
include("Perceptron_with_Potentials/CoupledBinaryPerceptrons.jl")
include("Perceptron_with_Potentials/CoupledBinaryPerceptronsNoTrace.jl")


include("Perceptron_with_Potentials/Potentials/TSGibbs/TSGibbs_RS.jl")
include("Perceptron_with_Potentials/Potentials/TSGibbs/TSGibbs_RSB.jl")
include("Perceptron_with_Potentials/Potentials/TSGibbs/TSGibbs_CoupledReplicas.jl")
include("Perceptron_with_Potentials/Potentials/TSGibbs/TSGibbs_CoupledReplicasNoTrace.jl")


include("Perceptron_with_Potentials/Potentials/TSPerceptron/TSPerceptron_RS.jl")
include("Perceptron_with_Potentials/Potentials/TSPerceptron/TSPerceptron_RSB.jl")
include("Perceptron_with_Potentials/Potentials/TSPerceptron/TSPerceptron_CoupledReplicas.jl")
include("Perceptron_with_Potentials/Potentials/TSPerceptron/TSPerceptron_CoupledReplicasNoTrace.jl")


include("Perceptron_with_Potentials/BinaryPerceptron1RSBDyntheta1withWick.jl")
include("Perceptron_with_Potentials/BinaryPerceptron1RSBDyntheta1NoWick.jl")


#RBMs with different units


# example model for docs
include("examples/HopfieldExample.jl")


include("FederatedHopfields/FederatedHopfield2.jl")
include("FederatedHopfields/FederatedHopfield2SameAlpha.jl")

include("Perceptron/SymmetricPerceptronTSQuenchedGibbsfixedR.jl")
include("Perceptron/SymmetricPerceptronTSQuenchedGibbsV1.jl")
include("Perceptron/SymmetricPerceptronTSQuenchedLinear.jl")

include("Perceptron/SymmetricPerceptronTSAnnealedGibbsMine.jl")
include("Perceptron/SymmetricPerceptronTSQuenchedGibbs.jl")

include("FederatedHopfields/v2_Corrupted_y2_RSB1_gianluca.jl")
include("FederatedHopfields/Corrupted.jl")
include("FederatedHopfields/HopfieldTS.jl")
include("FederatedHopfields/CoupledHopfieldTS.jl")
include("FederatedHopfields/Corrupted_TCpoint.jl")