/**
 * 
 */
package net.maizegenetics.analysis.monetdb;

import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Monetdb used LittleEndian format while java writes to BigEndian
 * This class writes the binary files to accepted monetdb format.
 * 
 * @author lcj34
 *
 */
public class LittleEndianDataOutputStream extends FilterOutputStream{

    public LittleEndianDataOutputStream(OutputStream out) {
        super(out);
    }
    
    public void writeByte(byte value) throws IOException{
        // Should be no translation needed
        out.write(value);
    }
    
    public void writeShort(short value) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(2).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putShort(value);
        out.write(buffer.array());
    }

    public void writeInt(int value) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putInt(value);
        out.write(buffer.array());
    }

    public void writeFloat(float value) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putFloat(value);
        out.write(buffer.array());
    }
    
    public void writeLong(long value) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(8).order(ByteOrder.LITTLE_ENDIAN);
        buffer.putLong(value);
        out.write(buffer.array());
    }
    
    public void writeChars(byte[] data) throws IOException {
        ByteBuffer buffer = ByteBuffer.allocate(data.length).order(ByteOrder.LITTLE_ENDIAN);
        buffer.put(data);
        out.write(buffer.array());
    }
//    public static void main(String... args) throws Exception {
//        short unsignedByte = 12;
//        int unsignedShort = 123;
//        long unsignedInt = 1234;
//
//        NativeDataOutputStream out =
//                new NativeDataOutputStream(
//                        new BufferedOutputStream(
//                                new FileOutputStream("writejava.bin", false)));
//
//        try {
//            out.write     ((byte)  unsignedByte );
//            out.writeShort((short) unsignedShort);
//            out.writeInt  ((int)   unsignedInt  );
//            out.flush();
//        } finally {
//            out.close();
//        }
//    }


}
