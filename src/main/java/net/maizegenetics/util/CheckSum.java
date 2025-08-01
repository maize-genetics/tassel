package net.maizegenetics.util;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.InputStream;
import java.security.MessageDigest;

public class CheckSum {

    private static final Logger myLogger = LogManager.getLogger(CheckSum.class);

    private CheckSum() {
        // utility class
    }

    public static String getMD5Checksum(String filename) {
        return getChecksum(filename, "MD5");
    }

    /**
     * Allows user to specify the protocol, e.g. MD5, SHA-1, SHA-256
     *
     * @param filename filename to checksum
     * @param protocol protocol
     *
     * @return check sum
     */
    public static String getProtocolChecksum(String filename, String protocol) {
        return getChecksum(filename, protocol);
    }

    /**
     * Allows user to specify the protocol, e.g. MD5, SHA-1, SHA-256
     *
     * @param filename filename to checksum
     * @param protocol protocol
     *
     * @return check sum
     */
    public static String getChecksum(String filename, String protocol) {

        try {
            InputStream inputStream = Utils.getInputStream(filename);
            MessageDigest digester = MessageDigest.getInstance(protocol);

            byte[] buffer = new byte[8192];
            int numOfBytesRead;
            while ((numOfBytesRead = inputStream.read(buffer)) > 0) {
                digester.update(buffer, 0, numOfBytesRead);
            }
            byte[] hashValue = digester.digest();
            return convertBytesToHex(hashValue);
        } catch (Exception ex) {
            myLogger.error(ex.getMessage());
        }

        return null;

    }

    private static String convertBytesToHex(byte[] bytes) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < bytes.length; i++) {
            builder.append(String.format("%02x", bytes[i] & 0xff));
        }
        return builder.toString();
    }

    /**
     * Allows user to specify the protocol, e.g. MD5, SHA-1, SHA-256
     *
     * @param str string to checksum
     * @param protocol protocol
     *
     * @return check sum
     */
    public static String getChecksumForString(String str, String protocol) {

        // from https://www.mkyong.com/java/java-md5-hashing-example/
        try {
            MessageDigest md = MessageDigest.getInstance(protocol);
            md.update(str.getBytes());
            byte byteData[] = md.digest();
            return convertBytesToHex(byteData);
        } catch (Exception e) {
            myLogger.error("getChecksumForString: problem getting checksum: " + e.getMessage());
            throw new IllegalStateException("CheckSum: getChecksumForString: error: " + e.getMessage());
        }

    }

}
