/** Register the four Nightingale components used by the offline structure viewer. */
import Linegraph from "@nightingale-elements/nightingale-linegraph-track";
import Navigation from "@nightingale-elements/nightingale-navigation";
import Sequence from "@nightingale-elements/nightingale-sequence";
import Track from "@nightingale-elements/nightingale-track";

for (const [name, component] of Object.entries({
  "linegraph-track": Linegraph, navigation: Navigation, sequence: Sequence, track: Track,
})) {
  const tag = `nightingale-${name}`;
  if (!customElements.get(tag)) customElements.define(tag, component);
}
