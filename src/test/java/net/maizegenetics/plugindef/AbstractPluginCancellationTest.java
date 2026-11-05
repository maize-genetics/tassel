package net.maizegenetics.plugindef;

import org.junit.Test;

import javax.swing.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * The progress panel's Cancel button interrupts the thread running a plugin, and the
 * plugin stops at its next progress report.
 */
public class AbstractPluginCancellationTest {

    private static class LoopingPlugin extends AbstractPlugin {

        private final CountDownLatch myStarted = new CountDownLatch(1);
        private final AtomicBoolean myFinished = new AtomicBoolean(false);

        LoopingPlugin() {
            super(null, false);
        }

        @Override
        public DataSet processData(DataSet input) {
            myStarted.countDown();
            long deadline = System.nanoTime() + TimeUnit.SECONDS.toNanos(30);
            while (System.nanoTime() < deadline) {
                fireProgress(50);
            }
            myFinished.set(true);
            return new DataSet(new Datum("result", "ok", null), this);
        }

        @Override
        public ImageIcon getIcon() {
            return null;
        }

        @Override
        public String getButtonName() {
            return "Looping Plugin";
        }

        @Override
        public String getToolTipText() {
            return "Looping Plugin";
        }

    }

    private static class RecordingListener implements PluginListener {

        private final AtomicBoolean myReturned = new AtomicBoolean(false);
        private final AtomicBoolean myFailed = new AtomicBoolean(false);

        @Override
        public void dataSetReturned(PluginEvent event) {
            myReturned.set(true);
        }

        @Override
        public void progress(PluginEvent event) {
        }

        @Override
        public void pluginFailed(Plugin plugin, Exception cause) {
            myFailed.set(true);
        }

    }

    @Test
    public void testInterruptStopsPluginAtNextProgressReport() throws Exception {
        LoopingPlugin plugin = new LoopingPlugin();
        RecordingListener listener = new RecordingListener();
        plugin.addListener(listener);

        PluginEvent event = new PluginEvent(new DataSet(new Datum[0], null));
        ThreadedPluginListener thread = new ThreadedPluginListener(plugin, event);
        thread.start();
        assertTrue("plugin never started", plugin.myStarted.await(10, TimeUnit.SECONDS));

        thread.interrupt();
        thread.join(TimeUnit.SECONDS.toMillis(10));

        assertFalse("plugin thread is still running", thread.isAlive());
        assertFalse("plugin ran to completion", plugin.myFinished.get());
        assertTrue(plugin.wasCancelled());
        assertFalse("cancelled plugin returned a data set", listener.myReturned.get());
        assertFalse("cancellation was reported as a failure", listener.myFailed.get());
    }

    @Test
    public void testPerformFunctionReturnsNullWhenCancelled() {
        AbstractPlugin plugin = new LoopingPlugin() {
            @Override
            public DataSet processData(DataSet input) {
                Thread.currentThread().interrupt();
                return super.processData(input);
            }
        };

        try {
            assertNull(plugin.performFunction(null));
            assertTrue(plugin.wasCancelled());
        } finally {
            Thread.interrupted();
        }
    }

}
