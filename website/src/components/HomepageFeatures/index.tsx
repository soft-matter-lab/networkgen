import type {ReactNode} from 'react';
import clsx from 'clsx';
import Heading from '@theme/Heading';
import styles from './styles.module.css';

type FeatureItem = {
  title: string;
  Svg?: React.ComponentType<React.ComponentProps<'svg'>>;
  img?: string;
  description: React.ReactNode;
};

const FeatureList: FeatureItem[] = [
  {
    title: 'Customize network architecture',
    Svg: require('@site/static/img/network_arch_main.svg').default,
    description: (
      <>
        NetworkGen allows for easy customization of network architecture, 
        including control over geometry, strand connectivity, and typology.
      </>
    ),
  },
  {
    title: 'Control network statistics',
    img: require('@site/static/img/network_gauss_main2.png').default,
    description: (
      <>
        Precisely control network statistics, including strand length distribution,
        crosslink functionality, and more.
      </>
    ),
  },
{
  title: 'Ready for LAMMPS',
  img: '/networkgen/img/networkfracture.gif',
  description: (
    <>
      NetworkGen data and table files are designed to work directly with LAMMPS. Generate your network
      and get started with your simulations! 
    </>
  ),
},
];

function Feature({title, Svg, img, description}: FeatureItem) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        {Svg 
          ? <Svg className={styles.featureSvg} role="img" />
          : <img src={img} alt={title} className={styles.featureSvg} />
        }
      </div>
      <div className="text--center padding-horiz--md">
        <Heading as="h3">{title}</Heading>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures(): ReactNode {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
