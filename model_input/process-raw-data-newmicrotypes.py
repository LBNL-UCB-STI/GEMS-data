import os
import shutil

import numpy as np
import pandas as pd


class ModeFile:
    def __init__(self, modeName, regions, dataFolder):
        self.__dataFolder = dataFolder
        self.__modeName = modeName
        self.__defaultColumns = ["MicrotypeID"]
        self.__df = pd.DataFrame(
            np.zeros((len(regions), 0)),
            index=regions,
        )

    def to_csv(self):
        modesDir = os.path.join(self.__dataFolder, "modes")
        if not os.path.exists(modesDir):
            os.makedirs(modesDir)
            print(modesDir)
        self.__df.to_csv(os.path.join(modesDir, self.__modeName + ".csv"))

    def adoptMissingColumns(self, otherDf):
        missingColumns = [col for col in otherDf.columns if col not in self.__defaultColumns]
        currentDf = self.__df.copy()
        for col in missingColumns:
            if pd.api.types.is_numeric_dtype(otherDf[col].dtype):
                self.__df[col] = otherDf[col].mean()
                currentDf[col] = otherDf[col].mean()


def mixedTrafficNetwork(lengths, MFDparams):
    MFDparams = MFDparams.set_index(["geotype", "network_microtype"]).reindex(lengths.index)
    df = lengths.to_frame()
    df["Dedicated"] = False
    df["MFD"] = "rural"
    df.loc[df.index.get_level_values(1) == "Urban_5", "MFD"] = "modified-quadratic"
    df["Type"] = "Road"
    df["avgLinkLength"] = 50
    # a,avgLinkLength,b,capacityFlow,k_jam,maxInflowDensity,maxInflowPerMeterPerHour,smoothingFactor,vMax,waveSpeed
    df["a"] = MFDparams["a"].copy()
    df["b"] = MFDparams["b"].copy()
    df["k_jam"] = 0.2
    df["maxInflowDensity"] = 0.15
    df["maxInflowPerMeterPerHour"] = 10.0
    df["vMax"] = MFDparams["max_spd_mph"] / 3600.0 * 1609.34
    return df


def walkNetwork(lengths):
    df = lengths.to_frame() / 5.0
    df["Dedicated"] = True
    df["MFD"] = "fixed"
    df["Type"] = "Sidewalk"
    df["vMax"] = 1.35
    return df


def busNetwork(lengths):
    df = lengths.to_frame() * 0.0
    df["Dedicated"] = True
    df["MFD"] = "bottleneck"
    df["Type"] = "Road"
    df["capacityFlow"] = 0.4
    df["vMax"] = 17.0
    df["avgLinkLength"] = 50
    return df


def railNetwork(lengths):
    df = lengths.to_frame() / 10.0
    df["Dedicated"] = True
    df["MFD"] = "fixed"
    df["Type"] = "Rail"
    df["vMax"] = 20.0
    return df


def createNetworks(inputDir, outputDir, diameter):
    if os.path.exists(os.path.join(inputDir, "SubNetworks.csv")):
        tots = pd.read_csv(
            os.path.join(inputDir, "SubNetworks.csv"),
            index_col=["geotype", "network_microtype"],
        )
        tots["Length"] = tots["laneMiles"] * 1609.34
    else:
        tractIDs = pd.read_csv(
            os.path.join(inputDir, "Typology", "microtype_geotype_output_2020.csv")
        )
        tractVals = pd.read_csv(
            os.path.join(inputDir, "Network", "network_microtype_metrics_2.csv")
        )
        joined = pd.merge(tractVals, tractIDs, left_on="tract", right_on="GEOID")
        tots = joined.groupby(["geotype", "network_microtype"]).agg({"laneMiles": "sum"})
        tots["Length"] = tots["laneMiles"] * 1609.34
    microtypes = pd.Series(
        diameter / 1609.34, index=tots.index, name="DiameterInMiles"
    ).to_frame()
    microtypes.to_csv(os.path.join(outputDir, "Microtypes.csv"))
    MFDparams = pd.read_csv(os.path.join(inputDir, "median_mfd_parameters_3x5.csv"))
    # CONVERT TO METERS AND SECONDS
    MFDparams["a"] *= (1609.34**2) / 3600
    MFDparams["b"] /= 1609.34

    subnetworks = {
        "Auto-Bus-Bike": mixedTrafficNetwork(tots["Length"], MFDparams),
        "Bus": busNetwork(tots["Length"]),
        "Walk": walkNetwork(tots["Length"]),
        "Rail": railNetwork(tots["Length"]),
    }

    allSubnetworks = pd.concat(subnetworks, axis=0, names=["Mode"]).reset_index()
    allSubnetworks.index.rename("SubnetworkID", inplace=True)
    print("stop")
    modeToSubnetwork = (
        allSubnetworks.Mode.replace(
            "Auto-Bus-Bike", "Auto-Bus-Bike-Freight_combi-Freight_single"
        )
        .str.lower()
        .str.split("-")
        .explode()
        .reset_index()
    )
    return allSubnetworks, microtypes, modeToSubnetwork


def defineLaneLengthAdjustment(microtypes, tps):

    out = pd.DataFrame(index=microtypes.index, columns=tps.index, dtype="object")
    cols = out.columns.copy()
    out = out.reset_index()
    for hr in cols:
        if int(hr) < 6:
            out[hr] = out["geotype"] + "-" + out["network_microtype"] + "-" + "overnight"
        elif int(hr) < 10:
            out[hr] = out["geotype"] + "-" + out["network_microtype"] + "-" + "am"
        elif int(hr) < 17:
            out[hr] = out["geotype"] + "-" + out["network_microtype"] + "-" + "midday"
        elif int(hr) < 21:
            out[hr] = out["geotype"] + "-" + out["network_microtype"] + "-" + "pm"
        else:
            out[hr] = out["geotype"] + "-" + out["network_microtype"] + "-" + "overnight"
    out = (
        out.set_index(["geotype", "network_microtype"])
        .stack()
        .rename("TimeMicrotypeGroupID")
        .to_frame()
    )
    out["LaneLengthAdjustment"] = 1.0
    return out


geotype = None
inFolder = "raw-data-national_3by5"
outFolder = "input-data-national_3by5"
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

newDir = os.path.join(ROOT_DIR, "..", "data", outFolder)
# if os.path.exists(newDir):
#     shutil.rmtree(newDir)
# os.makedirs(newDir)

occupancy = pd.read_csv(
    os.path.join(ROOT_DIR, "..", "data", inFolder, "avg_rider_per_trip.csv")
)
occupancy.rename(columns={"weighted_pax_per_trip": "Value"}, inplace=True)

occupancy["Mode"] = "auto"
occupancy.to_csv(os.path.join(ROOT_DIR, "..", "data", newDir, "ModeOccupancy.csv"))

subnetworks, ms, mts = createNetworks(
    os.path.join(ROOT_DIR, "..", "data", inFolder), newDir, 100
)

timePeriods = pd.read_csv(
    os.path.join(ROOT_DIR, "..", "data", inFolder, "TimePeriods.csv"), index_col="TimePeriodID"
)

lla = defineLaneLengthAdjustment(ms, timePeriods)
lla.to_csv(os.path.join(ROOT_DIR, "..", "data", newDir, "LaneLengthAdjustment.csv"))
subnetworks.to_csv(os.path.join(ROOT_DIR, "..", "data", newDir, "SubNetworks.csv"))
mts.to_csv(os.path.join(ROOT_DIR, "..", "data", newDir, "ModeToSubNetwork.csv"))

coeffs = pd.read_csv(
    os.path.join(ROOT_DIR, "..", "data", inFolder, "mode_choice_coefficients_implc_v2.csv")
)
coeffs = (
    coeffs.loc[coeffs["trip_purpose_agg"] != "social"]
    .replace({"shopping_meals": "leisure"})
    .rename(
        columns={
            "trip_purpose_agg": "TripPurposeID",
            "populationgroupid": "PopulationGroupID",
            "BetaTravelTimeDistBin1Bin2": "BetaTravelTime",
            "BetaTravelTimeDistBin3Plus": "BetaTravelTime_Long",
            "BetaMonetaryCostDistBin1Bin2": "BetaMonetaryCost",
            "BetaMonetaryCostDistBin3Plus": "BetaMonetaryCost_Long",
            "mode": "Mode",
        }
    )
)
coeffs["BetaTravelTime_Long"] -= coeffs["BetaTravelTime"]
coeffs["BetaMonetaryCost_Long"] -= coeffs["BetaMonetaryCost"]

coeffs["BetaTravelTimeMixed"] = 0.0
coeffs.loc[coeffs["Mode"] == "bike", "BetaTravelTimeMixed"] = (
    coeffs.loc[coeffs["Mode"] == "bike", "BetaTravelTime"].values / 4.0
)
coeffs.loc[coeffs["Mode"] == "bike", "BetaTravelTime"] *= 0.75
coeffs.to_csv(os.path.join(ROOT_DIR, "..", "data", newDir, "PopulationGroups.csv"), index=False)


for copiedFile in [
    "DistanceBins",
    "DistanceDistribution",
    "FreightDemand",
    "OriginDestination",
    "Population",
    "RoadNetworkCosts",
    "PopulationGroups",
    "TimePeriods",
    "TripGeneration",
    "TripPurposes",
    "Externalities",
]:
    oldPath = os.path.join(ROOT_DIR, "..", "data", inFolder, copiedFile + ".csv")
    newPath = os.path.join(ROOT_DIR, "..", "data", outFolder, copiedFile + ".csv")
    shutil.copyfile(oldPath, newPath)


oldPath = os.path.join(ROOT_DIR, "..", "data", inFolder, "TransitionMatrix-100m.csv")
newPath = os.path.join(ROOT_DIR, "..", "data", outFolder, "TransitionMatrix.csv")
shutil.copyfile(oldPath, newPath)


for mode in ["auto", "bike", "bus", "rail", "walk"]:
    oldPath = os.path.join(
        ROOT_DIR, "..", "data", "input-data-geotype-A", "modes", mode + ".csv"
    )
    newPath = os.path.join(ROOT_DIR, "..", "data", outFolder, "modes", mode + ".csv")
    df = pd.read_csv(oldPath)
    modeDf = ModeFile(mode, ms.index, newDir)
    modeDf.adoptMissingColumns(df)
    modeDf.to_csv()

ma = {}
for g_id, m_id in ms.index:
    new_index = ms.index.set_names(["d_geotype", "d_network_microtype"])
    df = pd.DataFrame(
        np.zeros((len(ms), 4)),
        index=new_index,
        columns=["bus_rail_bike", "bus_norail_bike", "nobus_rail_bike", "bus_rail_nobike"],
    )
    df["bus_rail_bike"] = 0.85
    df["nobus_rail_bike"] = 0.05
    df["bus_norail_bike"] = 0.05
    df["bus_rail_nobike"] = 0.05
    ma[(g_id, m_id)] = df

ma = pd.concat(ma, names=["o_geotype", "o_network_microtype"]).stack()
ma.index.rename("TransitLayer", level=4, inplace=True)
ma.rename("Portion").to_frame().to_csv(
    os.path.join(ROOT_DIR, "..", "data", newDir, "ModeAvailability.csv")
)
