import type {ReactNode} from 'react';
import clsx from 'clsx';
import Link from '@docusaurus/Link';
import useDocusaurusContext from '@docusaurus/useDocusaurusContext';
import useBaseUrl from '@docusaurus/useBaseUrl';
import Layout from '@theme/Layout';
import HomepageFeatures from '@site/src/components/HomepageFeatures';
import Heading from '@theme/Heading';
import styles from './index.module.css';

function HomepageHeader() {
  const {siteConfig} = useDocusaurusContext();
  const gifUrl = useBaseUrl('/img/networkfracture.gif');
  return (
<header className={styles.heroBanner}>
  <div className={styles.heroBannerInner}>
    <div className={styles.heroLeft}>
      <p className={styles.heroLabel}>Soft Matter Lab · University of Colorado Boulder</p>
      <Heading as="h1" className={styles.heroTitle}>
        A general-purpose polymer network generator
      </Heading>
      <p className={styles.heroSubtitle}>
        Precisely configure network geometry, strand length distributions, and defect
        structures. Generate output files ready for LAMMPS.
      </p>
      <div className={styles.buttons}>
        <Link className="button button--primary button--lg" to="/docs/overview">
          Get started
        </Link>
        <Link className={clsx('button button--lg', styles.btnOutline)} to="/config-builder">
          Config builder →
        </Link>
      </div>
    </div>
    <div className={styles.heroRight}>
      <img
        src={gifUrl}
        alt="NetworkGen fracture simulation"
        className={styles.heroGif}
      />
    </div>
  </div>
</header>
  );
}

export default function Home(): ReactNode {
  const {siteConfig} = useDocusaurusContext();
  return (
    <Layout
      title="NetworkGen — Polymer Network Generator"
      description="Generate simulation-ready polymer networks with precise control over geometry, strand length distributions, and defect structures.">
      <HomepageHeader />
      <main>
        <HomepageFeatures />
      </main>
    </Layout>
  );
}
