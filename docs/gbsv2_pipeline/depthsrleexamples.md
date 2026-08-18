The taxa/tag distribution data is stored in the GBSv2 database in the tagtaxadistribution table.  For memory efficienty reasons, this data is stored as a "blob" in the SQLite database.

The taxa are gathered from the key file when the GBSSeqToTagDBPlugin is run.  They are sorted alphabetically, then inserted to the Taxon table in alphabetic order. The taxonID in the Taxon table identifies these taxon.  The tagtaxadistribution table contains the taxon ID in the depthsRLE "blob" as well as depths information.

There may be times a user would like to view this data.  The data is encoded/decoded via the java classes TaxaDistribution and  TaxaDistBuilder, found in folder net.maizegenetics.dna.tag in the TASSEL source code.  Examples of coding/decoding data for the taxa distribution are below.

## Encoding taxa distribution:

Taxa distribution is created using the TaxaDistBuilder to create a TaxaDistribution object.  To create a TaxaDistribution with 75000 taxa but no values:

```
#!java

   TaxaDistribution myTD = TaxaDistBuilder.create(75000);
```

To increment by 1 the taxa count for a particular taxa id (e.g. taxa with id 25) do:

```
#!java

   myTD.increment(25);
```

To increment by 1 the taxa count for a series of taxa defined in a list:

```
#!java

  int[] taxaScored = {0,5,10,1005,2500,500,6289,70000};
  TaxaDistribution myTD = TaxaDistBuilder.create(75000);
  for (int index : taxaScored)
  {
     myTD.increment(index);
  }

```

To create a TaxaDistribution with room for 100 taxa and add to it a taxa list with specific depths:

```
#!java

  int[] taxaWithTags={1,5,8,9,10};
  int[] depthTags={1,1,2,313,2};
  TaxaDistribution myTD=TaxaDistBuilder.create(100,taxaWithTags,depthTags);
```
The encode TaxaDepth() method of the TaxaDistribution class is used to encode the TaxaDistribution as a "blob" of bytes for storing in the GBSv2 SQLite database table "tagtaxadistribution".

```
#!java

  byte[] encodeDepths=myTD.encodeTaxaDepth();
```

## Decoding taxa distribution from the GBSv2 database

Taxadistribution is stored in the tagtaxadistirubtion table of the GBSv2 database. In the source code, all methods to access the database are coded in class net.maizegenetics.dan.tag.TagDataSQLite.java.  New queries to the database should be coded in this method.  In addition to defined database access queries, this method keeps a tagTagIDMap for quicker access to a tag's id, as well as a list of taxa.

To get the taxa distribution for a particular tag, these steps may be followed:

```
#!java

// Get the tagID for the tag from the tagTagIDMap
int myTagid = tagTagIDMap.get(tag);  // "tag" is user supplied Tag object

// select only the depthsRLE from the tagtaxadistribution table
ResultSet rs = connection.createStatement().executeQuery("select depthsRLE from tagtaxadistribution where tagid=myTagid");

//Create the TaxaDistribution from the data returned
TaxaDistribution  myTD = TaxaDistBuilder.create(rs.getBytes(1));
```
Here is an example of creating a map of all tags with their depths for a particular taxon.  This example pulls both the tagID and the depths for a particular taxon.  It makes use of both data stored in the database and of the tagTagIDMap and taxaList maintained in the TagDataSQLite.java class.

```
#!java

    public Map<Tag,Integer> getTagsAndDepthForTaxon(Taxon taxon) {
        ImmutableMap.Builder<Tag,Integer> tdMapBuilder = new ImmutableMap.Builder<>();
        int taxonIndex=myTaxaList.indexOf(taxon);
        try {
            ResultSet rs=connection.createStatement().executeQuery("select * from tagtaxadistribution");
            while(rs.next()) {
                int tagDepth = TaxaDistBuilder.create(rs.getBytes("depthsRLE")).depths()[taxonIndex];
                if(tagDepth > 0) {
                    Tag addTag = tagTagIDMap.inverse().get(rs.getInt("tagid"));
                    tdMapBuilder.put(addTag,tagDepth);
                }
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        return tdMapBuilder.build();
    }
```
Please see TASSEL source code file net.maizegenetics.dna.tag.TagDataSQLite.java for additional examples of accessing data in the GBSv2 database.
