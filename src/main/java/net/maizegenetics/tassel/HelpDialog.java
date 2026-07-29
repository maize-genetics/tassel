/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package net.maizegenetics.tassel;

import javax.swing.*;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import java.awt.*;
import java.io.IOException;
import java.net.URL;

//Designed and programmed by Dr. Edward Buckler and his Bioinformatics Team: Peter Bradbury, Terry Casstevens, Chunguang Du, Dallas Kroon, Jack Liu, David Remington, Jeff Thornsberry, and Zhiwu Zhang.
public class HelpDialog extends JDialog {

    private JEditorPane htmlPane;

    public HelpDialog(Frame frame) {
        super(frame, "TASSEL Help", false);

        //Create the HTML viewing pane.
        htmlPane = new JEditorPane();
        htmlPane.setEditable(false);
        initHelp();
        JScrollPane htmlView = new JScrollPane(htmlPane);
        Dimension minimumSize = new Dimension(400, 400);
        htmlView.setMinimumSize(minimumSize);
        htmlView.setPreferredSize(new Dimension(600, 400));

        //Add the split pane to this frame.
        getContentPane().add(htmlView, BorderLayout.CENTER);

        pack();
    }

    private void initHelp() {
        String s = "Home.html";
        try {
            displayURL(s);
        } catch (Exception e) {
            System.err.println("Couldn't create net.maizegenetics.help URL: " + s);
        }
    }

    private void displayURL(String url) {
        try {
            htmlPane.setPage(HelpDialog.class.getResource(url));
            htmlPane.addHyperlinkListener(new HyperlinkListener() {
                public void hyperlinkUpdate(HyperlinkEvent hyperlinkEvent) {
                    HyperlinkEvent.EventType type = hyperlinkEvent.getEventType();
                    final URL url = hyperlinkEvent.getURL();
                    if (type == HyperlinkEvent.EventType.ENTERED) {
                        System.out.println("URL: " + url);
                    } else if (type == HyperlinkEvent.EventType.ACTIVATED) {
                        System.out.println("Activated");

                        //do some thing here
                        try {
                            Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + url);
                            //Runtime.getRuntime().exec("rundll32 url.dll,FileProtocolHandler " + "www.cnn.com");
                            //"C:\\Maize\\tassel\\src\\net\\maizegenetics\\help\\LoadSSR.gif");
                            //"C:\\Maize\\tassel\\src\\net\\maizegenetics\\help\\Overview.html");
                        } catch (Exception er) {
                            er.printStackTrace();
                        }

                    }
                }
            });

        } catch (IOException e) {
            System.err.println("Attempted to read a bad URL: " + url);
        }
    }
}
