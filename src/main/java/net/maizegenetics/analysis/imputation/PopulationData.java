package net.maizegenetics.analysis.imputation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Pattern;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.BitSet;

public class PopulationData {
	public String name;
	public String chromosome;
	public ArrayList<String> members;
	public BitSet snpIndex;
	public String parent1;
	public String parent2;
	public double contribution1;
	public double contribution2;
	public int Fgen;
	public double inbredCoef = -1;
	public GenotypeTable original;
	public GenotypeTable imputed;
	public byte[] alleleA;
	public byte[] alleleC;
	public boolean useSubpopulations;
	ArrayList<TaxaList> subpopulationGroups;
	ArrayList<BitSet> subpoplationSiteIndex;

	/**
	 * @param popFilename	the name of the file containing pedigree information for a group of populations
	 * @return	a HashMap of family names (keys) and associated PopulationData objects (values) containing information for the families in the pedigree file
	 */
	public static ArrayList<PopulationData> readPedigreeFile(String popFilename) {
		Pattern tab = Pattern.compile("\t");
		LinkedHashMap<String, PopulationData> familyMap = new LinkedHashMap<String, PopulationData>();
		String input = "";
		try {
			BufferedReader br = new BufferedReader(new FileReader(popFilename));
			br.readLine();
			while ((input = br.readLine()) != null) {
				String[] info = tab.split(input);
				if (info[0].length() > 0 && !info[0].equalsIgnoreCase("NA")) {
					PopulationData family = familyMap.get(info[0]);
					if (family == null) {
						family = new PopulationData ();
						family.name = info[0];
						family.members = new ArrayList<String>();
						family.members.add(info[1]);
						family.parent1 = info[2];
						family.parent2 = info[3];
						family.contribution1 = Double.parseDouble(info[4]);
						family.contribution2 = Double.parseDouble(info[5]);
						try {
							family.inbredCoef = Double.parseDouble(info[6]);
						} catch (Exception e) {}
						
						familyMap.put(info[0], family);
					}
					else family.members.add(info[1]);
				}
			}
			br.close();
		} catch (IOException e) {
			System.out.println("at: " + input);
			e.printStackTrace();
			System.exit(-1);
		} catch (NumberFormatException nfe) {
			System.out.println("at: " + input);
			nfe.printStackTrace();
			System.exit(-1);
		}
		return new ArrayList<PopulationData>(familyMap.values());
	}
	
        /**
         * @param popFilename   the name of the file containing pedigree information for a group of populations
         * @param includeParents        should the parents be included as family members?
         * @return      a HashMap of family names (keys) and associated PopulationData objects (values) containing information for the families in the pedigree file
         */
        public static ArrayList<PopulationData> readPedigreeFile(String popFilename, String chrname, boolean includeParents) {
                Pattern tab = Pattern.compile("\t");
                LinkedHashMap<String, PopulationData> familyMap = new LinkedHashMap<String, PopulationData>();
                String input = "";
                try {
                        BufferedReader br = new BufferedReader(new FileReader(popFilename));
                        br.readLine();
                        while ((input = br.readLine()) != null) {
                                String[] info = tab.split(input);
                                if (info[0].length() > 0 && !info[0].equalsIgnoreCase("NA")) {
                                        PopulationData family = familyMap.get(info[0]);
                                        if (family == null) {
                                                family = new PopulationData ();
                                                family.name = info[0];
                                                family.chromosome = chrname;
                                                family.parent1 = info[2];
                                                family.parent2 = info[3];
                                                family.members = new ArrayList<String>();
                                                if (includeParents) {
                                                    family.members.add(family.parent1);
                                                    family.members.add(family.parent2);
                                                }
                                                family.members.add(info[1]);
                                                family.contribution1 = Double.parseDouble(info[4]);
                                                family.contribution2 = Double.parseDouble(info[5]);
                                                try {
                                                        family.inbredCoef = Double.parseDouble(info[6]);
                                                } catch (Exception e) {}
                                                
                                                familyMap.put(info[0], family);
                                        }
                                        else family.members.add(info[1]);
                                }
                        }
                        br.close();
                } catch (IOException e) {
                        System.out.println("at: " + input);
                        e.printStackTrace();
                        System.exit(-1);
                } catch (NumberFormatException nfe) {
                        System.out.println("at: " + input);
                        nfe.printStackTrace();
                        System.exit(-1);
                }
                return new ArrayList<PopulationData>(familyMap.values());
        }
}
