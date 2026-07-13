/*
 *  GeneratedGuiBoilerplate
 */
package net.maizegenetics.plugindef;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Marks auto-generated plugin accessors and GUI-only hook methods (for example
 * {@code getIcon()}, {@code getButtonName()}, {@code getToolTipText()}) so that
 * code-coverage tooling can exclude this boilerplate from coverage reports.
 *
 * <p>The retention is {@link RetentionPolicy#RUNTIME} because Kover's IntelliJ
 * coverage engine only honors {@code annotatedBy} filters for annotations with
 * BINARY or RUNTIME retention.</p>
 *
 * @author Terry Casstevens
 */
@Retention(RetentionPolicy.RUNTIME)
@Target({ElementType.METHOD, ElementType.TYPE})
public @interface GeneratedGuiBoilerplate {
}
