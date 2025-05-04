

#include "collectdensity.h"
#include "collectenergy.h"
#include "collectdomaingrowth.h"
#include "collectpcf.h"
#include "collectsf.h"
#include "collectdgpcf.h"
#include "collectresidual.h"
#include "collectcenterofmass.h"

#define COLLECT_SCALAR \
CollectDensity, CollectEnergy, CollectDomainGrowth, \
CollectDGPCF, CollectPCF, CollectSF, CollectResidual, CollectCenterOfMass
#define COLLECT_COMPLEX CollectDensity, CollectEnergy, CollectCenterOfMass

