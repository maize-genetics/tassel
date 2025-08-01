package net.maizegenetics.analysis.gbs.v2;

import net.maizegenetics.analysis.gbs.Barcode;
import java.util.*;

/**
 * This is an implementation of a trie (prefix tree) in Java.
 * Supports opeations like searching a string, searching a prefix, searching by prefix etc.
 * @author Janu Verma
 * jv367@cornell.edu
 **/
public class BarcodeTrie{
    private TrieNode root;
    private Map<String, Barcode> barcodeInformation;

    /** Constructor */
    public BarcodeTrie() {
        root = new TrieNode();
        barcodeInformation = new HashMap<>();
    }

    /** Adds a Barcode to the trie
     * @param barcode
     */
    public void addBarcode(Barcode barcode){
        // Store both barcode and initial cut site
        String[] barcodeWOverhang = barcode.getBarWOverHang();
        for (String word: barcodeWOverhang) {
            root.addWord(word.toUpperCase());
            String bcode = word;
            barcodeInformation.put(bcode, barcode);
        }
    }
    
    /** Adds a Barcode to the trie
     * @param barcode
     */
    public void addBarcodeNoOverhang(Barcode barcode){
        // Store only barcode, NOT initial cut site
        String barcodeAlone = barcode.getBarcodeString();
        root.addWord(barcodeAlone.toUpperCase());
        String bcode = barcodeAlone;
        barcodeInformation.put(bcode, barcode);
 
    }

    /**
     * Add a collection of barcodes to the trie.
     * @param barcodes
     */
    public void addAllBarcodes(Collection<Barcode> barcodes){
        for (Barcode b: barcodes)
            addBarcode(b);
    }

    /**
     * Add a collection of barcodes with OUT overhangs to the trie.
     * @param barcodes
     */
    public void addAllBarcodesNoOverhang(Collection<Barcode> barcodes){
        for (Barcode b: barcodes)
            addBarcodeNoOverhang(b);
    }

    /**
     * checks if the String is in the trie.
     * @param s
     * @return true if the string is in the trie
     */
    public boolean contains(String s){
        TrieNode currentNode = root;
        for (int i = 0; i < s.length(); i++){
            char c = s.charAt(i);
            if (currentNode.containsKey(c))
                currentNode = currentNode.getNode(c);
            else
                return false;
        }
        return true;
    }



    /**
     * get the words in the trie with the given
     * prefix
     * @param prefix
     * @return a List contaning String objects containing the words
     * in the Trie with the given prefix.
     */
    public List getWords(String prefix){
        // Find the node which represents the last letter of the prefix.
        TrieNode lastNode = root;
        for (int i = 0; i < prefix.length(); i++){
            lastNode = lastNode.getNode(prefix.charAt(i));

            // If no node matches, then no words exits, return empty list
            if (lastNode == null) return new ArrayList();
        }
        // Return the words which eminate from the last node
        return lastNode.getWords();
    }

    /**
     * Find the longest prefix of a string
     * @param input
     */
    public Barcode longestPrefix(String input){
        String result = "";
        if(input==null) {
            System.out.println("stop");
        }
        int length = input.length();
        TrieNode crawl = root;
        int level, prevMatch = 0;
        for (level = 0; level < length-1; level++){
            char ch = input.charAt(level);
            if(ch<'A' || ch>'T') {
                ch=Character.toUpperCase(ch);
                if(ch<'A' || ch>'T') return null;
            }
            TrieNode child = crawl.getNode(ch);  //Get the Node reprsenting the character.
            if (crawl.containsKey(ch)){
                result += ch;
                crawl = child;
                if (crawl.isWord)
                    prevMatch = level + 1;
            }
            else break;
        }
        if (!crawl.isWord) result = result.substring(0,prevMatch);
        else result = result;
        return barcodeInformation.get(result);
    }


    public static void main(String args[]){
    }

    class TrieNode{
        public TrieNode parent;
        public TrieNode[] children;
        public boolean isLeaf; // Quick way to check if any children exist
        public boolean isWord; // does this node represent teh last character
        public char character; //character the node represents


        /**
         * Constructor for top level root node.
         */
        public TrieNode()
        {
            children = new TrieNode[26];
            isLeaf = true;
            isWord = false;
        }


        /**
         * Constructor for the child node.
         */
        public TrieNode(char character){
            this();
            this.character = character;
        }



        /**
         * Adds a word to this node. This method is called recursively and
         * adds child nodes for each successive letter in the word,
         * therefore recursive calls will be made with partial words.
         * @param word - the word to add
         */
        protected void addWord(String word){
            isLeaf = false;
            int charPos = word.charAt(0) - 'A';
            if (children[charPos] == null){
                children[charPos] = new TrieNode(word.charAt(0));
                children[charPos].parent = this;
            }
            if (word.length() > 1){
                children[charPos].addWord(word.substring(1));
            }
            else{
                children[charPos].isWord = true;
            }

        }




        /**
         * Return the child TrieNode representing rthe given char,
         * or null if no node exists.
         * @param c
         * @return TrieNode
         */
        protected TrieNode getNode(char c){
            return children[c-'A'];
        }



        /**
         * checks if the given character is a children.
         * @param c
         * @return true if c is a children
         */
        public boolean containsKey(char c){
            List followers = new ArrayList();
            for (TrieNode x : children) {
                if (x != null) {
                    char y = x.character;
                    followers.add(y);
                }
            }
            return (followers.contains(c));
        }



        /**
         * Returns a List of String objects which are lower in the
         * hierarchy than this node.
         * @return List of words
         */
        protected List getWords() {
            //Create a list to return.
            List list = new ArrayList();

            // If this node represents a word, add it.
            if (isWord) {
                list.add(toString());
            }
            // if any children
            if (!isLeaf) {
                //Add any words belonging to any children
                for (int i = 0; i < children.length; i++) {
                    if (children[i] != null) {
                        list.addAll(children[i].getWords());
                    }
                }
            }
            return list;
        }

        /**
         * Gets the string that this node represents.
         * e.g.g if this node represents the charcter t, whose parent
         * represents the character a, whose parent represents the charcter
         * c, then the string would be cat.
         * @return String
         */

        public String toString(){
            if (parent == null){
                return "";
            }
            else{
                return parent.toString() + new String(new char[]{character});
            }
        }
    }
}

