class Exports:
    """
    This class implement many methods used to export the Network object in different format.
    In particular the exports format will be methods-oriented (MaBoSS, Ginsim, cobrexa and so on...).
    To start with, the user can export the Network in SIF and Bnet format.
    In the future many more versatile methods will be implemented (SBML) and annotations will be included for each
    interaction, including the DOI of the relative reference and annotations from each database
    """
    def __init__(self, network):
        net = network.copy()
        net.convert_edgelist_into_genesymbol()
        self.nodes = net.nodes
        self.interactions = net.edges
        return

    def export_bnet(self, file_name="logic_model.bnet"):
        """
        Function to export the network in bnet format
        """
        with open(file_name, "w") as f:
            f.write("# model in BoolNet format\n")
            f.write("# the header targets, factors is mandatory to be importable in the R package BoolNet\n")
            f.write("\n")
            f.write(
                "targets, factors\n")  # this is the standard label that I found in the biolqm github, is it ok?
            for entry in self.nodes.values:
                node = entry[0]
                formula_on = self.interactions.loc[self.interactions["target"] == node].query("Effect == 'stimulation'")["source"].to_list()
                formula_off = self.interactions.loc[self.interactions["target"] == node].query("Effect == 'inhibition'")["source"].to_list()
                formula_undefined = self.interactions.loc[self.interactions["target"] == node].query("Effect == 'inhibition'")["source"].to_list()
                formula = formula_on + formula_off
                commons = list(set(formula_on).intersection(set(formula_off)))
                if formula_undefined:
                    print("The network has some undefined interaction that will be ignored. WARNING: this can result in"
                          "a disconnected model! ")
                    print(formula_undefined)
                # print(shared)
                for common in commons:
                    print("Two possible opposite interactions found for: ", common, " and ", node)
                    formula_off.remove(common) #here I should insert something to create an ensamble of model
                f.write(node + ",")
                offset = 16 - len(node)  # nice offset so the visualization is understandable
                f.write(" " * offset)
                if not formula:
                    f.write(" ( ")
                    f.write(node)
                    f.write(" ) ")
                    f.write("\n")
                if formula_on:
                    f.write(" ( ")
                    f.write(" | ".join(formula_on))  # writing the first parenthesis with all the positive interactions
                    f.write(" ) ")
                    if not formula_off:
                        f.write("\n")
                if formula_on != [] and formula_off != []:
                    f.write(" & ")
                    f.write(" !( ")
                    f.write(" | ".join(formula_off))  # writing the first parenthesis with all the positive and negative interactions
                    f.write(" ) ")
                    f.write("\n")
                if formula_on == [] and formula_off != []:
                    f.write(" !( ")
                    f.write(" | ".join(formula_off))  # writing the first parenthesis with all the negative interactions
                    f.write(" ) ")
                    f.write("\n")
        f.close  # good to go
        return

    def export_sif(self):
        """
        Function to export the network in SIF format
        """
        return
